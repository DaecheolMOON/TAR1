import argparse
import subprocess
import sys
import shutil
import tempfile
import os
import glob

def check_dependencies(tools):
    """필요한 커맨드라인 도구가 설치되어 있는지 확인합니다."""
    print("--- Checking for required tools... ---")
    all_found = True
    for tool in tools:
        if not shutil.which(tool):
            print(f"Error: Required tool '{tool}' is not installed or not in your PATH.", file=sys.stderr)
            all_found = False
    if not all_found:
        print("Please install the missing tools and try again.", file=sys.stderr)
        sys.exit(1)
    print("All required tools found.\n")

def process_tbl_to_bed_data(tbl_content):
    """
    grep, awk, sort 로직을 파이썬으로 직접 처리합니다.
    nhmmer의 .tbl 출력 내용을 받아 BED 형식 데이터 리스트를 반환합니다.
    """
    bed_lines = []
    for line in tbl_content.strip().split('\n'):
        # grep -v "^#"
        if line.startswith('#'):
            continue
        
        fields = line.split()
        if len(fields) < 10:
            continue
            
        # awk 'BEGIN{OFS="\t"} { a=$9; b=$10; start=(a<b?a:b)-1; end=(a<b?b:a); print $1, start, end }'
        target_name = fields[0]
        pos1 = int(fields[8])
        pos2 = int(fields[9])
        
        start = min(pos1, pos2) - 1
        end = max(pos1, pos2)
        
        # BED 포맷은 start가 0-based, end가 1-based 입니다.
        # 생성된 start, end 값은 이미 이 규칙을 따릅니다.
        bed_lines.append((target_name, start, end))

    # sort -k1,1 -k2,2n
    bed_lines.sort(key=lambda x: (x[0], x[1]))
    
    return [f"{name}\t{start}\t{end}" for name, start, end in bed_lines]

def run_command(command, step_name):
    """주어진 명령어를 실행하고 오류를 처리합니다."""
    print(f"--- Running: {step_name} ---")
    print(f"Command: {' '.join(command)}")
    try:
        result = subprocess.run(
            command, 
            check=True, 
            capture_output=True, 
            text=True
        )
        print(f"--- {step_name} successful. ---\n")
        return result
    except subprocess.CalledProcessError as e:
        print(f"\nError occurred during: {step_name}", file=sys.stderr)
        print(f"Command failed with exit code {e.returncode}", file=sys.stderr)
        print(f"Command: {' '.join(e.cmd)}", file=sys.stderr)
        print("\n--- STDOUT ---", file=sys.stderr)
        print(e.stdout, file=sys.stderr)
        print("\n--- STDERR ---", file=sys.stderr)
        print(e.stderr, file=sys.stderr)
        # 개별 파일 처리 실패 시 전체 스크립트를 중단하지 않고 다음 파일로 넘어갈 수 있도록
        # sys.exit(1) 대신 False를 반환합니다.
        return False
    except FileNotFoundError:
        print(f"Error: Command '{command[0]}' not found. Is it installed and in your PATH?", file=sys.stderr)
        sys.exit(1)

def run_pipeline(hmm_file, fasta_file, output_fasta_file, cpu, incE):
    """단일 FASTA 파일에 대한 전체 파이프라인을 실행합니다."""
    print(f"Processing file: {os.path.basename(fasta_file)}")
    
    # 임시 디렉토리 생성하여 중간 파일들을 관리
    with tempfile.TemporaryDirectory() as temp_dir:
        # 중간 파일 경로 설정
        temp_tbl_path = os.path.join(temp_dir, "hits.tbl")
        temp_bed_path = os.path.join(temp_dir, "blocks.bed")

        # 1. nhmmer 실행
        nhmmer_cmd = [
            'nhmmer',
            '--cpu', str(cpu),
            '--tblout', temp_tbl_path,
            '--incE', str(incE),
            hmm_file,
            fasta_file
        ]
        if not run_command(nhmmer_cmd, "nhmmer search"):
            return False # 실패 시 이 파일 처리를 중단

        # 2. tbl 파일 읽기 및 BED 데이터로 변환
        print("--- Running: Process TBL to BED ---")
        try:
            with open(temp_tbl_path, 'r') as f:
                tbl_content = f.read()
            
            bed_data = process_tbl_to_bed_data(tbl_content)
            
            with open(temp_bed_path, 'w') as f:
                f.write('\n'.join(bed_data))
            print("--- TBL to BED conversion successful. ---\n")
        except Exception as e:
            print(f"Error processing TBL file: {e}", file=sys.stderr)
            return False

        # 3. bedtools getfasta 실행
        # 최종 출력 디렉토리가 존재하는지 확인하고 생성
        os.makedirs(os.path.dirname(output_fasta_file), exist_ok=True)
        
        bedtools_cmd = [
            'bedtools', 'getfasta',
            '-fi', fasta_file,
            '-bed', temp_bed_path,
            # '-name', # getfasta는 기본적으로 BED의 4번째 필드를 이름으로 사용합니다. 
                      # 우리 BED 파일에는 3개 필드만 있으므로 이 옵션은 혼란을 줄 수 있습니다.
                      # 대신 좌표 기반으로 fasta 헤더를 생성하도록 합니다.
            '-fo', output_fasta_file
        ]
        if not run_command(bedtools_cmd, "bedtools getfasta"):
            return False

    print(f"✅ Successfully processed and saved to: {output_fasta_file}\n")
    return True

def main():
    """스크립트의 메인 로직을 실행합니다."""
    parser = argparse.ArgumentParser(
        description="A Python script to run a nhmmer -> bedtools pipeline on multiple fasta files in a directory.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '--hmm', 
        required=True, 
        help="Path to the input profile HMM file (e.g., TAR1.hmm)"
    )
    parser.add_argument(
        '--fasta', 
        required=True, 
        help="Path to the input directory containing reference FASTA files."
    )
    parser.add_argument(
        '--output', 
        required=True, 
        help="Path to the base directory for the output files."
    )
    parser.add_argument(
        '--cpu', 
        type=int, 
        default=6, 
        help="Number of CPU cores to use for nhmmer (default: 6)"
    )
    parser.add_argument(
        '--incE', 
        type=float, 
        default=1e-5, 
        help="Inclusion E-value threshold for nhmmer (default: 1e-5)"
    )
    args = parser.parse_args()

    # 0. 의존성 확인
    check_dependencies(['nhmmer', 'bedtools'])

    # 입력 FASTA 디렉토리 확인
    if not os.path.isdir(args.fasta):
        print(f"Error: Input FASTA path is not a directory: {args.fasta}", file=sys.stderr)
        sys.exit(1)

    # 출력 기본 디렉토리 생성
    os.makedirs(args.output, exist_ok=True)

    # FASTA 파일 목록 가져오기
    fasta_files = glob.glob(os.path.join(args.fasta, "*.fasta"))
    if not fasta_files:
        print(f"No .fasta files found in directory: {args.fasta}", file=sys.stderr)
        sys.exit(0)

    print(f"Found {len(fasta_files)} FASTA files to process.\n")
    
    success_count = 0
    fail_count = 0

    # 각 FASTA 파일에 대해 파이프라인 실행
    for fasta_path in fasta_files:
        filename = os.path.basename(fasta_path)
        
        # 파일 이름에서 샘플 ID 추출
        # e.g., "lengthresult_SRR26842280_with_quality.fasta" -> "SRR26842280"
        if filename.startswith("lengthresult_") and filename.endswith("_with_quality.fasta"):
            sample_id = filename.replace("lengthresult_", "").replace("_with_quality.fasta", "")
        else:
            # 예상치 못한 파일 이름 형식일 경우, 확장자만 제거
            sample_id = os.path.splitext(filename)[0]
            print(f"Warning: Filename '{filename}' does not match expected format. Using '{sample_id}' as sample ID.")

        # 출력 폴더 및 파일 경로 설정
        sample_output_dir = os.path.join(args.output, sample_id)
        final_output_path = os.path.join(sample_output_dir, "tar1_blocks.fa")
        
        # 파이프라인 실행
        if run_pipeline(args.hmm, fasta_path, final_output_path, args.cpu, args.incE):
            success_count += 1
        else:
            fail_count += 1
            print(f"Pipeline failed for {filename}.\n", file=sys.stderr)

    print("="*40)
    print("All tasks completed.")
    print(f"Total files processed: {len(fasta_files)}")
    print(f"✅ Successful: {success_count}")
    print(f"❌ Failed: {fail_count}")
    print("="*40)


if __name__ == "__main__":
    main()

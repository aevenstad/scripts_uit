# Import necessary modules
import os
import argparse
import sys
import subprocess
from Bio import SeqIO

# Define ANSI escape codes for text formatting
RED = '\033[91m'
GREEN = '\033[92m'
YELLOW = '\033[93m'
BOLD = '\033[1m'
RESET = '\033[0m'  # Reset text color to default


def check_existing_databases(fasta_files):
    existing_databases = []

    db_extensions = ['.nsd', '.nog', '.nhr', '.nin', '.nsq', '.nsi', '.psd', '.pog', '.phr', '.pin', '.psq', '.psi']

    for fasta_file in fasta_files:
        db_exists = any(os.path.exists(f"{fasta_file}{ext}") for ext in db_extensions)
        existing_databases.append((fasta_file, db_exists))

    return existing_databases


if __name__ == "__main__":
    # Create an argument parser
    parser = argparse.ArgumentParser(description="Make BLAST databases from FASTA files")
    
    # Add argument for processing all files in the current directory
    parser.add_argument("--all", action="store_true", help="Process all FASTA files in the current directory")
    
    # Add argument for specifying a specific file
    parser.add_argument("--in", dest="input_file", help="Specify a specific input FASTA file")
    
    args = parser.parse_args()

    # If argument --all is used
    if args.all:
        fasta_files = [file for file in os.listdir() if file.endswith(".fasta") or file.endswith(".fa")]
        existing_databases = check_existing_databases(fasta_files)
        
        print(f"\n##############################################################")
        print(f"################# {GREEN}Searching for fasta files{RESET} ##################")
        print(f"##############################################################\n")
        print(f"Looking for BLAST databases in {os.getcwd()}\n")

        # Print files with existing databases and files without databases
        print(f"{BOLD}Files with existing databases:{RESET}\n")
        for fasta_file, db_exists in existing_databases:
            if db_exists:
                print(f"- {fasta_file}")

        print(f"\n{BOLD}Files that need processing:{RESET}")

        files_needing_processing = []
        for fasta_file, db_exists in existing_databases:
            if not db_exists:
                print(f"- {fasta_file}")
                files_needing_processing.append(fasta_file)


        print(f"\n##############################################################")
        print(f"################## {GREEN}Building BLAST databases{RESET} ##################")
        print(f"##############################################################")
    

    # If argument --in is used
    elif args.input_file:
        # Process the specific input file
        fasta_files = [args.input_file]
        files_needing_processing = fasta_files
 

    
    # Process files without databases
    for fasta_file in files_needing_processing:
        print(f"\nProcessing: {RED}{fasta_file}{RESET}")
        dbtype = input("Enter database type: 'nucl' or 'prot':\n").strip()
        title = input("Enter a title for the database:\n").strip()                              

        command = f"makeblastdb -in {fasta_file} -dbtype '{dbtype}' -title '{title}' -parse_seqids"
        try:
            subprocess.run(command, check=True, shell=True)
        except subprocess.CalledProcessError as e:
            print(f"{YELLOW}Error!{RESET}")
            print(f"Command '{command}' failed with return code {e.returncode}\nRun command: 'conda activate blast'")
            print("Exiting")
            sys.exit(1)
        except FileNotFoundError:
            print(f"{YELLOW}Error!{RESET}")
            print(f"Command '{command}' not found")
            sys.exit(1)

if args.all:
    print(f"{GREEN}All databases created.{RESET}")
if args.input_file:
    print(f"{GREEN}Database created for {fasta_file}{RESET}")

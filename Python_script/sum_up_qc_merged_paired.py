# Author: Ugur Cabuk
#!/bin/python3
import os
import sys
def sum_up_num_reads(file1_path, file2_path, output_path):
    # Read data from the first file
    data1 = {}
    with open(file1_path, 'r') as file1:
        next(file1)  # skip the header
        for line in file1:
            fields = line.strip().split('\t')
            name = fields[0]
            num_reads = float(fields[4])
            data1[name] = num_reads
    
    with open(file2_path, 'r') as file2: #second file
        with open(output_path, 'w') as output_file:
            header = next(file2).strip()
            output_file.write(header + '\n') 
            for line in file2:
                fields = line.strip().split('\t')
                name = fields[0]
                num_reads = float(fields[4])
                if name in data1:
                    sum_num_reads = data1[name] + num_reads
                    fields[4] = str(sum_num_reads)
                output_file.write('\t'.join(fields) + '\n')


def process_directory(input_directory, output_directory):
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    for file_name in os.listdir(input_directory):
        if file_name.endswith("_merged"):
            merged_file_path = os.path.join(input_directory, file_name, "quant.sf")
            paired_file_path = os.path.join(input_directory, file_name.replace("_merged", "_paired"), "quant.sf")
            output_file_name = file_name.replace("_merged", "_merged_paired.quant.sf")
            output_file_path = os.path.join(output_directory, output_file_name)
            sum_up_num_reads(merged_file_path, paired_file_path, output_file_path)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_directory_path> <output_directory_path>")
        sys.exit(1)

    input_directory_path = sys.argv[1]
    output_directory_path = sys.argv[2]
    process_directory(input_directory_path, output_directory_path)


import csv

# Replace the paths with your own
txt_read_path = r"./temp.txt"
csv_write_path = r"./temp.csv"

with open(txt_read_path, "r") as txt_file:
    # Read all lines of the text file
    lines = txt_file.readlines()
    
    # Remove extra whitespaces and seperate by comma
    lines = [",".join(line.split()) for line in lines]

    # Write to a new CSV file
    with open(csv_write_path, "w") as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        for line in lines:
            writer.writerow(line.split(","))
#---------------------------------------------------------------------------------------------
# Parse the TkLayout cabling map files, to convert them to the format
# expected by this unpacker code.
#---------------------------------------------------------------------------------------------

import os

# Parse files for + & - ends of Tracker, to relate detid to PS or 2S module type.

def process_files_posneg(input_file_pos, input_file_neg):
    with open(input_file_pos, 'r') as f1, open(input_file_neg, 'r') as f2:
      # Concatenate files, skipping 3 lines of comments at top of each.
       combined_content = f1.readlines()[3:] + f2.readlines()[3:]

    mapDetType = {}
    for line in combined_content:
        parts = line.split(',')
        if len(parts) > 1: # Skip empty lines
          detid = int(parts[0])  
          type = parts[13].strip()
          if type == "2S":
            detType = "2S"
          elif "PS" in type:
            detType = "PS"
          else:
            print("An error has occurred"); exit(1)
          mapDetType[detid] = detType  # 2S or PS
        
    return mapDetType

# parse main file, to produce required cabling map output.
  
def process_file_main(input_file_main, input_file_pos, input_file_neg, output_file):

    # The Pos & Neg files are only needed to determine if each module is PS or 2S.
    mapDetType = process_files_posneg(input_file_pos, input_file_neg)
    
    # Read the main input file and collect the rows
    rows = []
    with open(input_file_main, 'r') as infile:
        # Skip first 3 lines of documentation
        next(infile); next(infile); next(infile)
      
        for line in infile:
            # Split the line by commas to get the three integers
            numbers = line.split(',')
            if len(numbers) > 1: # Skip empty lines
              det_id = int(numbers[0])
              dtc_id = int(numbers[1])
              dtc_ch = int(numbers[2])

              # Determine module type (PS or 2S)
              det_type = mapDetType[det_id]
            
              # Append a tuple of the four integers to the rows list
              rows.append((det_id, dtc_id, dtc_ch, det_type))
    
    # Sort the rows based on the third integer, then by the second integer if third integers are identical
    rows.sort(key=lambda x: (x[2], x[1]))
    
    # Write the sorted rows to the output file
    with open(output_filename, 'w') as outfile:
        for aa, bb, cc, dd  in rows:
            
            # Write the result to the output file
            outfile.write(f'{aa}  {cc}  {bb}  {dd}\n')

#---------- Main ----------

# IAN: TO DO -- persuade tkLayout authors to put this info into a single CSV file.
input_filename_main = 'TkLayoutFiles/CMSSWCablingMapOuter.csv' 
input_filename_pos= 'TkLayoutFiles/ModulesToDTCsPosOuter.csv' 
input_filename_neg = 'TkLayoutFiles/ModulesToDTCsNegOuter.csv' 
output_filename = './tkLayoutCabling_parsed.txt'

if os.path.exists(output_filename):
  os.remove(output_filename)

process_file_main(input_filename_main, input_filename_pos, input_filename_neg,
                             output_filename)

print("File created:", output_filename)

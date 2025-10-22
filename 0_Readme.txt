This release contains files necessary to run Sekhmet and examples of results.
Necessary files are :
   - A stellar model
   - "input.par" : list of all input parameters and their value (with indication of the reference value for each)
   - Sekhmet.x (compiled version of "Sekhmet_Release-1.f90")
   - Amon-Ra_Release1.py

Example files are :
   - "Sekhmet_output.csv" : result file for the reference scenario
   - "Sekhmet_output_a1.50.csv" : result file for a scenario with an initial semi-major axis of 1.5 au

How to run Sekhmet :
   - Define the initial parameters in the "input.par" file
   - Run Sekhmet.x
   - It produces output files "Sekhmet_output[...].csv" and a log file "Sekhmet.log".
   - Run Amon-Ra with the output file(s) that you want to consider to prin graphs of "a", "e" or "dot(a)" (edit the .py file to change the printed parameter, default is "a")

Contact:
Antoine Thuillier
antoine.thuillier@uliege.be
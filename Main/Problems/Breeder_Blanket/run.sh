cd ../..
trelis mesher.py "json_input='3-MOT/Main/Problems/Breeder_Blanket/Parameters/MOT_parameters_breeder_blankets.json'" 
trelis mesher.py  json_input=''

cd 3-MOT/Main/Problems/Breeder_Blanket/Source_terms/

fenicsproject run 

python3 source_term_writer.py #check hard coded json file

cd 3-MOT/Main


python3 3-MOT.py #check hard coded json file
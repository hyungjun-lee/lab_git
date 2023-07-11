# o2-analysis-jetfinderQA --aod-file @input_data.txt | 
o2-analysis-je-jet-finder-qa-hj - --configuration json://configuration.json -b @input_data.txt | 
# o2-jet-tutorial-skeleton - --aod-file @input_data.txt | 
# o2-analysistutorial-jet-analysis -b --aod-file @input_data.txt |
#o2-analysis-lf-lambdakzeroanalysis-mc  --configuration json://configuration.json -b   |
#o2-analysistutorial-jet-analysis - --configuration json://configuration.json -b |
#o2-analysis-lf-lambdakzerobuilder - --configuration json://configuration.json  -b  |
#o2-analysis-event-selection - --configuration json://configuration.json  -b  |
#o2-analysis-multiplicity-table  --configuration json://configuration.json -b   --aod-file @input_data.txt |
#o2-analysis-pid-tpc  --configuration json://configuration.json  -b  |
o2-analysis-event-selection - --configuration json://configuration.json  -b |
o2-analysis-timestamp - --configuration json://configuration.json  -b |
o2-analysis-track-propagation - --configuration json://configuration.json  -b |
# o2-analysis-collision-converter - --configuration json://configuration.json  -b  |
o2-analysis-alice3-trackselection - --configuration json://configuration.json  -b |
o2-analysis-je-jet-finder - --configuration json://configuration.json  -b |
o2-analysis-je-jet-matching  - --configuration json://configuration.json  -b |
o2-analysis-track-to-collision-associator - --configuration json://configuration.json  -b 
#o2-analysis-mm-track-propagation  -b --configuration json://configuration.json   --aod-file @input_data.txt |
#o2-analysis-trackselection - --configuration json://configuration.json  -b  |
#o2-analysis-centrality-table - --configuration json://configuration.json  -b  | 
#o2-analysis-hf-candidate-creator-2prong - --configuration json://configuration.json  -b  | 
#o2-analysis-je-jet-finder-hf - --configuration json://configuration.json  -b  | 
#o2-analysis-je-jet-matching-hf - --configuration json://configuration.json  -b |
#o2-analysis-hf-refit-pv-dummy - --configuration json://configuration.json  -b  |
#o2-analysis-hf-candidate-creator-3prong - --configuration json://configuration.json  -b  
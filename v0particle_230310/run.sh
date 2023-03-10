################### No options
#o2-analysis-mm-dndeta-hi --configuration json://configuration.json -b |

#o2-analysis-event-selection - --configuration json://configuration.json  -b |
#o2-analysis-track-propagation - --configuration json://configuration.json  -b |
#o2-analysis-collision-converter - --configuration json://configuration.json  -b | 
#o2-analysis-timestamp  --configuration json://configuration.json  -b

################### V0particle Wo cent
o2-analysis-mm-dndeta-hi --configuration json://configuration.json -b |
o2-analysis-mm-track-propagation  -b --configuration json://configuration.json   --aod-file @input_data.txt |

o2-analysis-event-selection - --configuration json://configuration.json  -b |
o2-analysis-track-propagation - --configuration json://configuration.json  -b |
o2-analysis-collision-converter - --configuration json://configuration.json  -b | 
o2-analysis-timestamp --configuration json://configuration.json  -b |
o2-analysis-trackselection - --configuration json://configuration.json  -b |
o2-analysis-pid-tpc  --configuration json://configuration.json  -b  |
o2-analysis-pid-tof  --configuration json://configuration.json  -b  |
o2-analysis-multiplicity-table  --configuration json://configuration.json -b   --aod-file @input_data.txt |
o2-analysis-pid-tof-base - --configuration json://configuration.json  -b  |
o2-analysis-pid-tpc-base - --configuration json://configuration.json  -b  |

o2-analysis-lf-lambdakzerobuilder - --configuration json://configuration.json  -b  

################### V0particle W cent
# o2-analysis-mm-track-propagation  -b --configuration json://configuration.json   --aod-file @input_data.txt |

# o2-analysis-event-selection - --configuration json://configuration.json  -b |
# o2-analysis-track-propagation - --configuration json://configuration.json  -b |
# o2-analysis-collision-converter - --configuration json://configuration.json  -b | 
# o2-analysis-timestamp --configuration json://configuration.json  -b |
# o2-analysis-trackselection - --configuration json://configuration.json  -b |
# o2-analysis-pid-tpc  --configuration json://configuration.json  -b  |
# o2-analysis-pid-tof  --configuration json://configuration.json  -b  |
# o2-analysis-pid-tpc-base - --configuration json://configuration.json  -b  | 
# o2-analysis-multiplicity-table  --configuration json://configuration.json -b   --aod-file @input_data.txt |
# o2-analysis-pid-tof-base - --configuration json://configuration.json  -b  |
# o2-analysis-pid-tpc-full - --configuration json://configuration.json  -b  |
# o2-analysis-centrality-table - --configuration json://configuration.json  -b |

# o2-analysis-lf-lambdakzerobuilder - --configuration json://configuration.json  -b  |
# o2-analysis-mm-dndeta-hi --configuration json://configuration.json -b 

################### D0particle tagging /// debuging...
#o2-analysis-mm-dndeta-hi --configuration json://configuration.json -b |
#o2-analysis-mm-track-propagation  -b --configuration json://configuration.json   --aod-file @input_data.txt |

#o2-analysis-event-selection - --configuration json://configuration.json  -b |
#o2-analysis-track-propagation - --configuration json://configuration.json  -b |
#o2-analysis-collision-converter - --configuration json://configuration.json  -b | 
#o2-analysis-timestamp --configuration json://configuration.json  -b |
#o2-analysis-trackselection - --configuration json://configuration.json  -b |
#o2-analysis-multiplicity-table  --configuration json://configuration.json -b   --aod-file @input_data.txt |
#o2-analysis-pid-tpc  --configuration json://configuration.json  -b  |
#o2-analysis-pid-tpc-full - --configuration json://configuration.json  -b  |
#o2-analysis-pid-tof-full - --configuration json://configuration.json  -b  |
#o2-analysis-pid-tof-base - --configuration json://configuration.json  -b  |

#o2-analysis-lf-lambdakzerobuilder - --configuration json://configuration.json  -b |

#o2-analysis-hf-track-index-skim-creator  - --configuration json://configuration.json  -b |
#o2-analysis-hf-candidate-creator-2prong - --configuration json://configuration.json  -b |
#o2-analysis-hf-candidate-selector-d0 - --configuration json://configuration.json  -b  


#o2-analysis-lf-lambdakzeroanalysis-mc  --configuration json://configuration.json -b   |
#o2-analysis-hf-task-bplus --configuration json://configuration.json -b   |
#o2-analysis-lf-lambdakzerobuilder - --configuration json://configuration.json  -b  |
#o2-analysis-event-selection - --configuration json://configuration.json  -b  |
#o2-analysis-multiplicity-table  --configuration json://configuration.json -b   --aod-file @input_data.txt |
#o2-analysis-pid-tpc  --configuration json://configuration.json  -b  |
#o2-analysis-timestamp  --configuration json://configuration.json  -b |
#o2-analysis-track-propagation - --configuration json://configuration.json  -b |
#o2-analysis-collision-converter - --configuration json://configuration.json  -b  |
#o2-analysis-mm-track-propagation  -b --configuration json://configuration.json   --aod-file @input_data.txt |
#o2-analysis-trackselection - --configuration json://configuration.json  -b |
#o2-analysis-pid-tpc-full - --configuration json://configuration.json  -b  |
#o2-analysis-pid-tof-full - --configuration json://configuration.json  -b  |
#o2-analysis-pid-tof-base - --configuration json://configuration.json  -b  |
#o2-analysis-hf-candidate-creator-2prong - --configuration json://configuration.json  -b  |
#o2-analysis-hf-track-index-skim-creator  - --configuration json://configuration.json  -b  |
#o2-analysis-hf-candidate-creator-bplus  - --configuration json://configuration.json  -b  |
#o2-analysis-hf-tree-creator-bplus-to-d0-pi - --configuration json://configuration.json  -b  |
#o2-analysis-mc-converter - --configuration json://configuration.json  -b |
#o2-analysis-hf-candidate-selector-d0 - --configuration json://configuration.json  -b 
#o2-analysis-centrality-table - --configuration json://configuration.json  -b   
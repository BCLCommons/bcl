This folder contains *example* descriptor files grouped by application

See aa/AADescriptorsHelp.txt or qsar/QSARDescriptorsHelp.txt complete descriptions of each descriptor 

Help can also be obtained for specific descriptors using GenerateDataset, e.g. for small molecule / qsar descriptors
  bcl.exe GenerateDataset -source SdfFile -feature_labels '3DA(help)'
or for proteins:
  bcl.exe GenerateDataset -source SequenceDirectory -feature_labels 'HelixStrandCoil(help)'

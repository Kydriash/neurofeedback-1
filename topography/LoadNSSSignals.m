function NSS = LoadNSSSignals(sNSSFilename);


x = xmlread(sNSSFilename);
[pathstr,name,ext] = fileparts(sNSSFilename);
cd(pathstr);
k=1;
vsChanNames = struct;
sChanNames  = x.getElementsByTagName('vsChanNames').item(0);

childNode = sChanNames.getFirstChild;
while ~isempty(childNode)
  %Filter out text, comments, and processing instructions.
  if childNode.getNodeType == childNode.ELEMENT_NODE
      %Assume that each element has a single org.w3c.dom.Text child
      childText = char(childNode.getFirstChild.getData);
      switch char(childNode.getTagName)
          case 'nChannels',
              vsChanNames.nChannels = childText;
          case 's' 
              vsChanNames.Name{k} = childText;
              k=k+1;
      end
  end
  childNode = childNode.getNextSibling;
end

NSS.vsChanNames = vsChanNames;

k=1;
vSignals = struct;
Signals= x.getElementsByTagName('vSignals');
for i=0:Signals.getLength-1
  thisSignal = Signals.item(i);
  childNode = thisSignal.getFirstChild;
  while ~isempty(childNode)
      %Filter out text, comments, and processing instructions.
      if childNode.getNodeType == childNode.ELEMENT_NODE
          %Assume that each element has a single org.w3c.dom.Text child
          childText = char(childNode.getFirstChild.getData);
          switch char(childNode.getTagName)
              case 'nSignals',
                  vSignals.nSignals = childText;
              case 'DerivedSignal' 
                  vSignals.DerivedSignal{k}.sSignalName = char(childNode.getElementsByTagName('sSignalName').item(0).item(0).get('TextContent'));
                  vSignals.DerivedSignal{k}.nChannels = str2num(childNode.getElementsByTagName('nChannels').item(0).item(0).get('TextContent'));
                  vSignals.DerivedSignal{k}.nBeamspaceDim = str2num(childNode.getElementsByTagName('nBeamspaceDim').item(0).item(0).get('TextContent'));
                  SpFilterFileName = (childNode.getElementsByTagName('SpatialFilterMatrix').item(0).item(0).get('TextContent'));
                  vSignals.DerivedSignal{k}.SpatialFilterMatrix = ReadEMSEMatrixXML(SpFilterFileName);
                  vSignals.DerivedSignal{k}.fBandpassLowHz = str2num(childNode.getElementsByTagName('fBandpassLowHz').item(0).item(0).get('TextContent'));
                  vSignals.DerivedSignal{k}.fBandpassHighHz = str2num(childNode.getElementsByTagName('fBandpassHighHz').item(0).item(0).get('TextContent'));
                  vSignals.DerivedSignal{k}.bHilbertEnvelope = str2num(childNode.getElementsByTagName('bHilbertEnvelope').item(0).item(0).get('TextContent'));
                  vSignals.DerivedSignal{k}.nLxNorm = str2num(childNode.getElementsByTagName('nLxNorm').item(0).item(0).get('TextContent'));
                  vSignals.DerivedSignal{k}.sTransform = char(childNode.getElementsByTagName('sTransform').item(0).item(0).get('TextContent'));
                  k=k+1;
          end
      end
      childNode = childNode.getNextSibling;
  end
 end
NSS.vSignals = vSignals; 
NSS.fSampleRateHz  = str2num(x.getElementsByTagName('fSampleRateHz').item(0).item(0).get('TextContent'));


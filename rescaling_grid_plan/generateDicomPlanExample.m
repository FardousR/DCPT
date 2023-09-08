function generateDicomPlanExample()

% Program made by PRP for manipulation of proton dicom plans

% This part inputs an eight layer mouse plan and outputs a plan with
% selected energy layers each containing a 10 x 10 PBS field
% These plans were used for scintillator sheets tests and for measuring
% depth dose curves in December 2022.

protonEnergyPerLayer = [150 120 80];
nLayersOut = length(protonEnergyPerLayer);

MUperSpot = [19 12 4.6];

spotSpacing = 2.5;
xPositions = -50:spotSpacing:50;
yPositions = -50:spotSpacing:50;

xPlanNew = [];
yPlanNew = [];
sign = 1;
for i = 1:length(xPositions)
    xPlanNew = [xPlanNew xPositions(i)*ones(size(yPositions))];
    yPlanNew = [yPlanNew sign*yPositions];
    sign = -sign;
end

% Template plan used to house the new plan:
planDir = 'C:\Data\PhD MSc projekter\Mouse experiments DCPT\Mouse plans';
fileNameIn = [planDir '\Plan22.dcm'];

meta = dicominfo(fileNameIn);
I = dicomread(meta);
metaOutput = meta;

fieldIDs = fieldnames(meta.('IonBeamSequence'));

% clear cumMetersetWeightPlan refDoseRefSeq protonEnergy numOfScanSpotPositions scanSpotSize scanSpotMetersetWeights scanSpotPositionMapX scanSpotPositionMapY

for fieldIndex = 1:length(fieldIDs)
    FieldNumber(fieldIndex) = meta.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('BeamNumber');
    meterSetRate(fieldIndex) =  meta.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence').('Item_1').('MetersetRate');
    finalCumulativeMetersetWeight(fieldIndex) = meta.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('FinalCumulativeMetersetWeight');
    NumOfCP(fieldIndex) = meta.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('NumberOfControlPoints');
    RefPatientSetupNumber(fieldIndex) = meta.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('ReferencedPatientSetupNumber');
    Isocenter(fieldIndex,1:3) = meta.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence').('Item_1').('IsocenterPosition')(1:3);
    
    cpIDs(fieldIndex,1:NumOfCP(fieldIndex)) = fieldnames(meta.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence'));
end

for cpIndex = 1:NumOfCP(fieldIndex)
    cumMetersetWeightPlan(fieldIndex,cpIndex) = meta.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence').(char(cpIDs(fieldIndex,cpIndex))).('CumulativeMetersetWeight');
    if( isfield(meta.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence').(char(cpIDs(fieldIndex,cpIndex))), 'NominalBeamEnergy'))
        protonEnergy(fieldIndex,cpIndex) =  meta.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence').(char(cpIDs(fieldIndex,cpIndex))).('NominalBeamEnergy');   % MeV
    end
    numOfScanSpotPositions(fieldIndex,cpIndex) = meta.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence').(char(cpIDs(fieldIndex,cpIndex))).('NumberOfScanSpotPositions');
    scanSpotSize(fieldIndex,cpIndex,1:2) = meta.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence').(char(cpIDs(fieldIndex,cpIndex))).('ScanningSpotSize')(1:2);
    
    scanSpotMetersetWeights(fieldIndex,cpIndex,1:numOfScanSpotPositions(fieldIndex,cpIndex)) = double(meta.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence').(char(cpIDs(fieldIndex,cpIndex))).('ScanSpotMetersetWeights')(1:numOfScanSpotPositions(fieldIndex,cpIndex)));
    scanSpotPositionMapX(fieldIndex,cpIndex,1:numOfScanSpotPositions(fieldIndex,cpIndex)) = meta.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence').(char(cpIDs(fieldIndex,cpIndex))).('ScanSpotPositionMap')(1:2:2*numOfScanSpotPositions(fieldIndex,cpIndex));
    scanSpotPositionMapY(fieldIndex,cpIndex,1:numOfScanSpotPositions(fieldIndex,cpIndex)) = meta.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence').(char(cpIDs(fieldIndex,cpIndex))).('ScanSpotPositionMap')(2:2:2*numOfScanSpotPositions(fieldIndex,cpIndex));
end

for fieldIndex = 1:length(fieldIDs)
    if( isfield(meta.('FractionGroupSequence').('Item_1').('ReferencedBeamSequence').(char(fieldIDs(fieldIndex))), 'BeamMeterset'))
        refBeamNumber = meta.('FractionGroupSequence').('Item_1').('ReferencedBeamSequence').(char(fieldIDs(fieldIndex))).('ReferencedBeamNumber');
        MUforField(RefPatientSetupNumber==refBeamNumber) = meta.('FractionGroupSequence').('Item_1').('ReferencedBeamSequence').(char(fieldIDs(fieldIndex))).('BeamMeterset');
    end
end
MUperCumMeterSetWeight =  MUforField./finalCumulativeMetersetWeight;


nLayersIn = NumOfCP(fieldIndex)/2;  % 8 layers

% Maniputale spots in the first already existing eight layers:
for layerNum = 1:min(nLayersIn,nLayersOut)
    cpIndex = 2*layerNum - 1;
    xPlan = squeeze(scanSpotPositionMapX(fieldIndex,cpIndex,1:numOfScanSpotPositions(fieldIndex,cpIndex)));
    yPlan = squeeze(scanSpotPositionMapY(fieldIndex,cpIndex,1:numOfScanSpotPositions(fieldIndex,cpIndex)));
    muPlan = squeeze(MUperCumMeterSetWeight(fieldIndex)*scanSpotMetersetWeights(fieldIndex,cpIndex,1:numOfScanSpotPositions(fieldIndex,cpIndex)));
    
    metaOutput.('RTPlanLabel') = 'Scintil tests';
    metaOutput.IonBeamSequence.Item_1.TreatmentMachineName = 'TR4';
    metaOutput.('ApprovalStatus') = 'APPROVED';
    
    numSpotsNew = length(xPlanNew);
    spotPositionMapNew = zeros(2*numSpotsNew,1);
    spotPositionMapNew(1:2:2*numSpotsNew) = xPlanNew;
    spotPositionMapNew(2:2:2*numSpotsNew) = yPlanNew;
    
    muPlanNew = MUperSpot(layerNum)*ones(size(xPlanNew));
    
    
    metaOutput.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence').(char(cpIDs(fieldIndex,cpIndex))).('ScanSpotPositionMap') = spotPositionMapNew;
    metaOutput.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence').(char(cpIDs(fieldIndex,cpIndex))).('NumberOfScanSpotPositions') =  numSpotsNew;
    metaOutput.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence').(char(cpIDs(fieldIndex,cpIndex))).('ScanSpotMetersetWeights') = muPlanNew/MUperCumMeterSetWeight(fieldIndex);
    metaOutput.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence').(char(cpIDs(fieldIndex,cpIndex))).('NominalBeamEnergy') = protonEnergyPerLayer(layerNum);
    
    metaOutput.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence').(char(cpIDs(fieldIndex,cpIndex+1))).('ScanSpotPositionMap') = spotPositionMapNew;
    metaOutput.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence').(char(cpIDs(fieldIndex,cpIndex+1))).('NumberOfScanSpotPositions') =  numSpotsNew;
    metaOutput.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence').(char(cpIDs(fieldIndex,cpIndex+1))).('ScanSpotMetersetWeights') = zeros(size(muPlanNew));
    metaOutput.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence').(char(cpIDs(fieldIndex,cpIndex+1))).('NominalBeamEnergy') = protonEnergyPerLayer(layerNum);
    
    % Set table/couch position to your experiment to avoid tedious overrides:
    metaOutput.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence').('Item_1').('TableTopVerticalPosition') =  -140;
    metaOutput.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence').('Item_1').('TableTopLongitudinalPosition') = 14;
    metaOutput.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence').('Item_1').('TableTopLateralPosition') = -356;
    metaOutput.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence').('Item_1').('PatientSupportAngle') = 234;
end



if nLayersOut > nLayersIn   % Create new control points to host the remaining energy layers:
    
    for cpIndex = NumOfCP(fieldIndex)+1:2:2*nLayersOut
        layerNum = (cpIndex+1)/2;
        
        % Update cpIDs with the new control point
        newFieldID = strcat('Item_',num2str(cpIndex));
        cpIDs(fieldIndex,cpIndex) = cellstr(newFieldID);
        
        % Use previous control point as a template for the new cp:
        newFieldStructure = metaOutput.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence').(char(cpIDs(fieldIndex,cpIndex-1)));
        newFieldStructure.('ScanSpotPositionMap') = spotPositionMapNew;
        newFieldStructure.('NumberOfScanSpotPositions') = numSpotsNew;
        newFieldStructure.('ScanSpotMetersetWeights') = muPlanNew(layerNum,:)/MUperCumMeterSetWeight(fieldIndex);
        newFieldStructure.('NominalBeamEnergy') = protonEnergyPerLayer(layerNum);
        [metaOutput.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence').(newFieldID)] = newFieldStructure;
        
        
        % Update cpIDs with the new control point
        newFieldID = strcat('Item_',num2str(cpIndex+1));
        cpIDs(fieldIndex,cpIndex+1) = cellstr(newFieldID);
        
        % Use previous control point as a template for the new cp:
        newFieldStructure = metaOutput.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence').(char(cpIDs(fieldIndex,cpIndex-1)));
        newFieldStructure.('ScanSpotPositionMap') = spotPositionMapNew;
        newFieldStructure.('NumberOfScanSpotPositions') = numSpotsNew;
        newFieldStructure.('ScanSpotMetersetWeights') = zeros(size(muPlanNew(layerNum,:)));
        newFieldStructure.('NominalBeamEnergy') = protonEnergyPerLayer(layerNum);
        [metaOutput.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence').(newFieldID)] = newFieldStructure;
    end    
elseif nLayersOut < nLayersIn  % Remove excess energy layers:
    for cpIndex = 2*nLayersOut+1:NumOfCP(fieldIndex)
        metaOutput.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence') = rmfield(metaOutput.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence'), (char(cpIDs(fieldIndex,cpIndex))));
    end
end


%**********************************************
% Adjust other parts of the Dicom plan structure to the changes made:
% This adjustment could be put into a separate function later on, which
% makes sure that the plan is made self consistent after manipulation such
% as addition/removal of spots or layers, or changing of energies
%**********************************************
cpIDsForField = fieldnames(metaOutput.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence'));
numOfCPsForField = length(cpIDsForField);
clear cumMetersetWeight
currentMetersetWeight = 0;
% Run through all layers and spots to find cumulative meterset weight
for cpIndex = 1:2:numOfCPsForField  % Each layer (defined by a specific proton energy) has two control points
    cumMetersetWeight(cpIndex) = currentMetersetWeight;
    ScanSpotMetersetWeightsForCP =  metaOutput.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence').(char(cpIDsForField(cpIndex))).('ScanSpotMetersetWeights');
    currentMetersetWeight = currentMetersetWeight + sum(ScanSpotMetersetWeightsForCP);
    cumMetersetWeight(cpIndex+1) = currentMetersetWeight;
    metaOutput.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence').(char(cpIDsForField(cpIndex))).('NumberOfScanSpotPositions') = length(ScanSpotMetersetWeightsForCP);
    metaOutput.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence').(char(cpIDsForField(cpIndex+1))).('NumberOfScanSpotPositions') = length(ScanSpotMetersetWeightsForCP);
end

% Set parameters for layer level
refDoseRefSeqForField = cumMetersetWeight/cumMetersetWeight(end);
for cpIndex = 1:numOfCPsForField
    metaOutput.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence').(char(cpIDsForField(cpIndex))).('ControlPointIndex') = cpIndex-1;
    metaOutput.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence').(char(cpIDsForField(cpIndex))).('CumulativeMetersetWeight') = cumMetersetWeight(cpIndex);
    metaOutput.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('IonControlPointSequence').(char(cpIDsForField(cpIndex))).('ReferencedDoseReferenceSequence').('Item_1').('CumulativeDoseReferenceCoefficient') = refDoseRefSeqForField(cpIndex);
end

% Set parameters for field level:
metaOutput.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('FinalCumulativeMetersetWeight') = cumMetersetWeight(end);  % This is MU/Gy as specified for the field in Eclipse
metaOutput.('IonBeamSequence').(char(fieldIDs(fieldIndex))).('NumberOfControlPoints') = numOfCPsForField;

% With the following setting the MUs (or the normalization) should be correct in Eclipse if the manipulated plan is normalized to the same value/percentage as the original plan after import into Eclipse
metaOutput.('FractionGroupSequence').('Item_1').('ReferencedBeamSequence').(char(fieldIDs(fieldIndex))).('BeamMeterset') =  MUperCumMeterSetWeight(fieldIndex)*cumMetersetWeight(end);

fileNameOut = [fileNameIn(1:end-10) 'ScintillatorTest_10x10_' num2str(protonEnergyPerLayer(1))];
for layerIndex = 2:nLayersOut
    fileNameOut = [fileNameOut '-' num2str(protonEnergyPerLayer(layerIndex))];
end
fileNameOut = [fileNameOut 'MeV.dcm'];
dicomwrite(I,fileNameOut, metaOutput, 'CreateMode', 'copy');









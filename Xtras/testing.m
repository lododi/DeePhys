

%% Filter the Units
units_FR = zeros(1,length(BT_rg.Units));
for iUnit = 1:length(BT_rg.Units)
    units_FR(iUnit) = height(BT_rg.Units(iUnit).SpikeTimes)/BT_rg.Units(iUnit).MEArecording.RecordingInfo.Duration;
end

BT_rg.Units = BT_rg.Units(units_FR>0.5);

%% Filter the Units within the Recordings
for iRec = 1:length(BT_rg.Recordings)
    mask_FR = zeros(1, length(BT_rg.Recordings(iRec).Units));
    for iUnit = 1:length(BT_rg.Recordings(iRec).Units)
        mask_FR(iUnit) = height(BT_rg.Recordings(iRec).Units(iUnit).SpikeTimes)/BT_rg.Recordings(iRec).RecordingInfo.Duration;
    end
    BT_rg.Recordings(iRec).Units = BT_rg.Recordings(iRec).Units(mask_FR>0.5);
end

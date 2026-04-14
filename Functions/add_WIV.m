function add_WIV(RecordingGroup)
    for k = 1:length(RecordingGroup.Recordings)
    
        DIV = RecordingGroup.Recordings(k).Metadata.DIV;
    
        WIV = floor((DIV - 2) / 7) + 1;
    
        RecordingGroup.Recordings(k).Metadata.WIV = WIV;
        RecordingGroup.Recordings(k).saveObject();
    end
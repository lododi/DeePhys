function DeePUOI = phy2DeePUOI(phyUOI,obj)
nSpikes_UOI = spike_units(spike_units==phyUOI);
spikes_timing = {obj.Units(1,:).SpikeTimes};
mask = cellfun(@(v) numel(v) == nSpikes_UOI, spikes_timing);
DeePUOI = find(mask==1);


nSpikes_UOI = length(obj.Spikes.Units(obj.Spikes.Units==phyUOI));
spikes_timing = {obj.Units(1,:).SpikeTimes};
mask = cellfun(@(v) numel(v) == nSpikes_UOI(1), spikes_timing);
DeePUOI = find(mask==1);
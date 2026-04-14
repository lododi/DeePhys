function convertedPoint = UMAPcoordinateconverter(pattern_rg,point)
arguments
    pattern_rg
    point
end
    UMAPmin = 125.75;
    UMAPmax = 375.25;
    TrueMin = min(pattern_rg.DimensionalityReduction.Recording.UMAP.Reduction); 
    TrueMax = max(pattern_rg.DimensionalityReduction.Recording.UMAP.Reduction);
    UMAPRange = UMAPmax - UMAPmin;
    TrueRange = TrueMax - TrueMin;
    convertedPoint = [(point(1,1)-UMAPmin)/(UMAPRange/TrueRange(1,1))+TrueMin(1,1),(point(1,2)-UMAPmin)/(UMAPRange/TrueRange(1,2))+TrueMin(1,2)];
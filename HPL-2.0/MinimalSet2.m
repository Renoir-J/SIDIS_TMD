HPL[{-1, 0}, t_]:=-HPL[{-2}, t] + HPL[{-1}, t]*HPL[{0}, t]/;(Not[$HPLAutoProductExpand]||$HPLAutoConvertToKnownFunctions)&&$HPLAutoReduceToMinimalSet;
HPL[{0, 0}, t_]:=HPL[{0}, t]^2/2/;(Not[$HPLAutoProductExpand]||$HPLAutoConvertToKnownFunctions)&&$HPLAutoReduceToMinimalSet;
HPL[{1, -1}, t_]:=HPL[{-1}, t]*HPL[{1}, t] - HPL[{-1, 1}, t]/;(Not[$HPLAutoProductExpand]||$HPLAutoConvertToKnownFunctions)&&$HPLAutoReduceToMinimalSet;
HPL[{1, 0}, t_]:=HPL[{0}, t]*HPL[{1}, t] - HPL[{2}, t]/;(Not[$HPLAutoProductExpand]||$HPLAutoConvertToKnownFunctions)&&$HPLAutoReduceToMinimalSet;
HPL[{1, 1}, t_]:=HPL[{1}, t]^2/2/;(Not[$HPLAutoProductExpand]||$HPLAutoConvertToKnownFunctions)&&$HPLAutoReduceToMinimalSet;
HPL[{-1, -1}, t_]:=HPL[{-1}, t]^2/2/;(Not[$HPLAutoProductExpand]||$HPLAutoConvertToKnownFunctions)&&$HPLAutoReduceToMinimalSet;

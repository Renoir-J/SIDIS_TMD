HPL[{minus, 0}, t_]:=HPL[{0}, t]*HPL[{minus}, t] - HPL[{0, minus}, t]/;(Not[$HPLAutoProductExpand]||$HPLAutoConvertToKnownFunctions)&&$HPLAutoReduceToMinimalSet;
HPL[{0, 0}, t_]:=HPL[{0}, t]^2/2/;(Not[$HPLAutoProductExpand]||$HPLAutoConvertToKnownFunctions)&&$HPLAutoReduceToMinimalSet;
HPL[{plus, minus}, t_]:=HPL[{minus}, t]*HPL[{plus}, t] - HPL[{minus, plus}, t]/;(Not[$HPLAutoProductExpand]||$HPLAutoConvertToKnownFunctions)&&$HPLAutoReduceToMinimalSet;
HPL[{plus, 0}, t_]:=HPL[{0}, t]*HPL[{plus}, t] - HPL[{0, plus}, t]/;(Not[$HPLAutoProductExpand]||$HPLAutoConvertToKnownFunctions)&&$HPLAutoReduceToMinimalSet;
HPL[{plus, plus}, t_]:=HPL[{plus}, t]^2/2/;(Not[$HPLAutoProductExpand]||$HPLAutoConvertToKnownFunctions)&&$HPLAutoReduceToMinimalSet;
HPL[{minus, minus}, t_]:=HPL[{minus}, t]^2/2/;(Not[$HPLAutoProductExpand]||$HPLAutoConvertToKnownFunctions)&&$HPLAutoReduceToMinimalSet;

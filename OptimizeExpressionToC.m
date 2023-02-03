(* ::Package:: *)

BeginPackage[ "OptimizeExpressionToC`"]

OptimizeExpressionToC::usage = "Generates optimized version of expression in C";

ExtractTags::usage = "";
ClearTagText::usage = "";
TagExistQ::usage = "";
AppendToTag::usage = "";

Begin[ "Private`"]

OptimizeExpressionToC[expr_List] :=
	Module[ {optimizedExpr, mainExpr, n, m, defs, output},
		optimizedExpr = Experimental`OptimizeExpression[expr];
		If[ ToString@optimizedExpr[[1, 0]] == "Block",
		    {n = Length[optimizedExpr[[1, 1]]];
		     mainExpr = optimizedExpr[[1, 2, n + 1]];},
		    {n = 0;
		     mainExpr = Flatten@{optimizedExpr[[1]]};}
		];
		m = Length[mainExpr];
		
		defs  = 
		Table[ "double " <>
			       ToString@CForm@optimizedExpr[[1, 2, i, 1]]  <>
			       " = " <>
			       ToString@CForm@optimizedExpr[[1, 2, i, 2]] <>
     			       ";",
    		       {i, 1, n}];
		
		output = 

		       MapIndexed[ "out("<>StringJoin@Riffle[ToString/@(#2-1),""]<>") = " <> ToString@CForm@#1 <>";" &,
					 mainExpr, {ArrayDepth[mainExpr]}];
		
		StringReplace[Join[defs,output], "Compile_$" -> "t"]
	];

End[]


ExtractTags[string_String] :=
	Module[ {regex, extractTag},
		regex = RegularExpression["\\n([^\\n]*)// (\\[|</?)([^\\]\\n]*)(\\]|>)[^\\n]*\\n"];
		
		extractTag[bounds_] := Module[{substring, indentation, tagType, tag},
					      substring = StringTake[string, bounds];
					      {{indentation, tagType, tag}} = 
					      StringCases[substring, regex -> {"$1", "$2", "$3"}];
					      Association["indentation" -> indentation,
							  "tagType" -> tagType, 
							  "tag" -> tag , "start" -> bounds[[1]] + 1, 
							  "end" -> bounds[[2]]]
				       ];
		
		extractTag[#] & /@ StringPosition[string, regex]
	];

TagExistQ[string_String, tag_String] := AnyTrue[ExtractTags[string], #[["tag"]]==tag &];

ClearTagText[string_String, tag_String] :=
	Module[{tags, startTag, endTag, x},
	       tags = ExtractTags[string];
	       startTag = FirstCase[tags, x_ /; x[["tag"]] == tag  && x[["tagType"]] == "<"];
	       endTag = FirstCase[tags, x_ /; x[["tag"]] == tag  && x[["tagType"]] == "</"];
	       
	       StringTake[string, {1, startTag[["end"]]}] <> 
			 StringTake[string, {endTag[["start"]], StringLength[string]}]
	]

AppendToTag[string_String,  tag_String, linesToAppend_List] := 
	Module[{tags, startTag, indent, x},
	       tags = ExtractTags[string];
	       startTag = FirstCase[tags, x_ /; x[["tag"]] == tag  && x[["tagType"]] == "<"];
	       indent = startTag[["indentation"]];
	       
	       StringTake[string, {1, startTag[["end"]]}] <> 
			 StringRiffle[ linesToAppend, {indent, "\n" <> indent, "\n"}] <>
			 StringTake[ string, {startTag[["end"]] + 1, StringLength[string]}]
	];

AppendToTag[string_String,  tag_String, stringToAppend_String] :=
	AppendToTag[string, tag, StringSplit[stringToAppend,"\n"]];


EndPackage[]




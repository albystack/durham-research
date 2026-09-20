(* Optional legacy renderer for the binary tiling matrix. Open this file in
   Mathematica, adjust the two paths if necessary, and evaluate the script. *)

scriptDirectory = DirectoryName[$InputFileName];
matrixPath = FileNameJoin[
  {scriptDirectory, "..", "output", "random_uniform_n200_seed20260728",
   "tiling_matrix.txt"}
];
outputPath = FileNameJoin[
  {scriptDirectory, "..", "output", "random_uniform_n200_seed20260728",
   "tiling_mathematica.png"}
];

AztecPrinter[x0_] := Module[{pieces = {}, size, i, j, centre, colour},
  size = Length[x0];
  For[i = 1, i <= size, i++,
    For[j = 1, j <= size, j++,
      If[x0[[i, j]] == 1,
        centre = {j - i, size + 1 - (i + j)};
        colour = Which[
          OddQ[i] && OddQ[j], Green,
          OddQ[i] && EvenQ[j], Blue,
          EvenQ[i] && EvenQ[j], Yellow,
          True, Red
        ];
        pieces = Join[pieces, {
          colour,
          If[OddQ[i] == OddQ[j],
            Rectangle[centre - {2, 1}, centre + {2, 1}],
            Rectangle[centre - {1, 2}, centre + {1, 2}]
          ]
        }]
      ]
    ]
  ];
  Rotate[pieces, -45 Degree]
];

tiling = Import[matrixPath, "Table"];
graphic = Graphics[
  AztecPrinter[tiling],
  Background -> White,
  ImageSize -> 1800
];
Export[outputPath, graphic];
graphic

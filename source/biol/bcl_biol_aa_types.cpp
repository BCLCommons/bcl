// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "biol/bcl_biol_aa_types.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

    //! @brief construct all AATypes
    AATypes::AATypes() :
//                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          Properties                                       charge pK     pk    pk    pk     pk    pk     pk     pk    pk      pk   pk     pk    pK     pk     pk     pk    pk   SSEProb.      Free Energies           Transfer Free Energies                                                  Free Energies continued
//                                                                   3L    1L  Natural Parent                   Atoms                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   type of first SCA   Prev   Steric  Pol     Vol   Hydrop. Isoe.              EMBOSS DTA   Sol   Sil    Rod   Pat    Wik    Leh   Gri     Bje  ProM   BjeN  BjeC   ProN   ProC   CarN  CarC     H      S      H       S       C       WW      ESG     KD      EISEN   HW      GUY     JANIN   PM1D    PM3D    Core    Trans   Sol   CoreH   TransH  SolH    CoreS   TransS  SolC    CoreC   TransC  SolC  CorePore CoreMem BBHydro   ExtBl  ExtTy     SASA  Girth SCPolr, SCTPSA,  VdwSA,  HAcc HDon, Aromatic Vdw Radius of C-beta
      ALA( AddEnum( "ALANINE",               /*  0 */   AATypeData( "ALA", 'A', true,  "ALA", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().H,   GetAtomTypes().HA,  GetAtomTypes().HB1, GetAtomTypes().HB2, GetAtomTypes().HB3),                                                                                                                                                                                                                                                                                                         GetAtomTypes().CB,  0.085, 1.280, 0.050, 1.000,  0.310,  6.110,/*A*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 7.59, 3.55,  7.58,  3.75,  9.69, 2.34, 0.420, 0.230, -0.168,  0.097,  0.116,  0.500, -1.600,  1.800,  0.620,  0.500,  0.100, -0.300, -0.170, -0.150, -0.130,  0.075,  0.081, -0.207, -0.071, -0.136, -0.109,  0.215,  0.211, -0.044,  0.177,  0.153, -2.588,  2.327,  1.454, -0.259, -0.597, 209.020, 2.15,  5.361,  40.07,  63.35,    0,   0, false, 1.61))), //ALA
      ARG( AddEnum( "ARGININE",              /*  1 */   AATypeData( "ARG", 'R', true,  "ARG", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().CG,  GetAtomTypes().CD,  GetAtomTypes().NE,  GetAtomTypes().CZ,  GetAtomTypes().NH1,  GetAtomTypes().NH2,  GetAtomTypes().H,    GetAtomTypes().HA,   GetAtomTypes().HB2,  GetAtomTypes().HB3,  GetAtomTypes().HG2,  GetAtomTypes().HG3,   GetAtomTypes().HD2,  GetAtomTypes().HD3,   GetAtomTypes().HE,  GetAtomTypes().HH11, GetAtomTypes().HH12, GetAtomTypes().HH21, GetAtomTypes().HH22), GetAtomTypes().CB,  0.051, 2.340, 0.290, 6.130, -1.010, 10.740,/*R*/    1.0, 12.5, 12.0, 12.5, 12.0, 11.50, 11.2, 12.48, 12.40, 12.0, 12.00, 12.50, 7.50, 3.55, 11.50, 11.50,  9.04, 2.17, 0.360, 0.250,  0.049, -0.002, -0.043,  1.810, 12.300, -4.500, -2.530, -3.000,  1.910,  1.400,  0.370,  0.320,  0.444, -0.128, -0.147,  0.927, -0.144, -0.185, -0.020, -0.061, -0.002,  0.156, -0.027, -0.097, -0.624, -4.596, -4.274,  0.000,  0.000, 335.732, 7.36,  2.222,   0.00,  42.20,    0,   3, false, 1.84))), //ARG
      ASN( AddEnum( "ASPARAGINE",            /*  2 */   AATypeData( "ASN", 'N', true,  "ASN", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().CG,  GetAtomTypes().OD1, GetAtomTypes().ND2, GetAtomTypes().H,   GetAtomTypes().HA,   GetAtomTypes().HB2,  GetAtomTypes().HB3,  GetAtomTypes().HD21, GetAtomTypes().HD22),                                                                                                                                                                                                                    GetAtomTypes().CB,  0.043, 1.600, 0.130, 2.950, -0.600,  6.520,/*N*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 7.50, 3.55,  7.22,  3.64,  8.80, 2.02, 0.210, 0.220,  0.157,  0.081, -0.182,  0.850,  4.800, -3.500, -0.780, -0.200,  0.480,  0.500,  0.180,  0.220,  0.397, -0.143, -0.114,  0.759,  0.055, -0.048,  0.048, -0.220,  0.241,  0.113, -0.243, -0.175, -2.750, -2.167, -2.367,  0.122, -0.210, 259.845, 4.35, 12.577,  63.64, 150.80,    1,   1, false, 1.78))), //ASN
      ASP( AddEnum( "ASPARTIC_ACID",         /*  3 */   AATypeData( "ASP", 'D', true,  "ASP", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().CG,  GetAtomTypes().OD1, GetAtomTypes().OD2, GetAtomTypes().H,   GetAtomTypes().HA,   GetAtomTypes().HB2,  GetAtomTypes().HB3,  GetAtomTypes().HD2),                                                                                                                                                                                                                                          GetAtomTypes().CB,  0.058, 1.600, 0.110, 2.780, -0.770,  2.950,/*D*/   -1.0,  3.9,  4.4,  3.9,  4.0,  3.68,  4.2,  3.90,  3.86,  3.5,  4.05,  4.07, 7.50, 4.55,  3.57,  4.57,  9.60, 1.88, 0.250, 0.200,  0.169,  0.067, -0.180,  3.640,  9.200, -3.500, -0.900, -3.000,  0.780,  0.600,  0.370,  0.410,  0.574, -0.068, -0.237,  0.963,  0.227, -0.193,  0.237, -0.295,  0.074,  0.266, -0.176, -0.253, -4.159, -4.524, -2.153,  0.181, -0.366, 257.993, 4.00,  6.164,  43.09,  75.41,    2,   0, false, 1.77))), //ASP
      CYS( AddEnum( "CYSTEINE",              /*  4 */   AATypeData( "CYS", 'C', true,  "CYS", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().SG,  GetAtomTypes().H,   GetAtomTypes().HA,  GetAtomTypes().HB2, GetAtomTypes().HB3,  GetAtomTypes().HG),                                                                                                                                                                                                                                                                                     GetAtomTypes().CB,  0.013, 1.770, 0.130, 2.430,  1.540,  6.350,/*C*/   -1.0,  8.5,  8.5,  8.3,  9.0,  8.33,  0.0,  8.18,  8.33,  6.8,  9.00,  8.28, 7.50, 3.55,  8.00,  9.00,  8.18, 1.96, 0.170, 0.410, -0.030,  0.245, -0.148, -0.020, -2.000,  2.500,  0.290,  1.000, -1.420, -0.900, -0.060, -0.150,  0.040,  0.369, -0.250, -0.125,  0.240, -0.285,  1.381,  1.205, -0.159, -0.228,  0.281, -0.353, -1.989,  0.378,  4.301,  0.094,  0.013, 240.500, 3.71,  5.222,  38.80,  56.60,    0,   0, false, 1.74))), //CYS
      GLN( AddEnum( "GLUTAMINE",             /*  5 */   AATypeData( "GLN", 'Q', true,  "GLN", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().CG,  GetAtomTypes().CD,  GetAtomTypes().OE1, GetAtomTypes().NE2, GetAtomTypes().H,    GetAtomTypes().HA,   GetAtomTypes().HB2,  GetAtomTypes().HB3,  GetAtomTypes().HG2,  GetAtomTypes().HG3,  GetAtomTypes().HE21, GetAtomTypes().HE22),                                                                                                                                                     GetAtomTypes().CB,  0.038, 1.560, 0.180, 3.950, -0.220,  5.650,/*Q*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 7.50, 3.55,  6.73,  3.57,  9.13, 2.17, 0.360, 0.250,  0.028, -0.046,  0.020,  0.770,  4.100, -3.500, -0.850, -0.200,  0.950,  0.700,  0.260,  0.030,  0.381, -0.039, -0.201,  0.797, -0.012, -0.270, -0.023, -0.136, -0.082,  0.170,  0.109, -0.056, -1.367, -2.601, -3.214,  0.042, -0.295, 286.761, 5.50,  7.999,  43.09,  94.99,    1,   1, false, 1.81))), //GLN
      GLU( AddEnum( "GLUTAMIC_ACID",         /*  6 */   AATypeData( "GLU", 'E', true,  "GLU", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().CG,  GetAtomTypes().CD,  GetAtomTypes().OE1, GetAtomTypes().OE2, GetAtomTypes().H,    GetAtomTypes().HA,   GetAtomTypes().HB2,  GetAtomTypes().HB3,  GetAtomTypes().HG2,  GetAtomTypes().HG3,  GetAtomTypes().HE2),                                                                                                                                                                           GetAtomTypes().CB,  0.069, 1.560, 0.150, 3.780, -0.640,  3.090,/*E*/   -1.0,  4.1,  4.4,  4.3,  4.5,  4.25,  4.2,  4.07,  4.25,  4.2,  4.45,  4.45, 7.70, 4.75,  4.15,  4.75,  9.67, 2.19, 0.420, 0.210,  0.013,  0.080, -0.082,  3.630,  8.200, -3.500, -0.740, -3.000,  0.830,  0.700,  0.150,  0.300,  0.553,  0.071, -0.318,  0.905,  0.090, -0.418,  0.167, -0.024, -0.039,  0.262,  0.139, -0.271, -1.607, -4.323, -2.144,  0.009, -1.768, 285.025, 5.05,  7.196,  40.07,  82.94,    2,   0, false, 1.84))), //GLU
      GLY( AddEnum( "GLYCINE",               /*  7 */   AATypeData( "GLY", 'G', true,  "GLY", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().H,  GetAtomTypes().HA2, GetAtomTypes().HA3),                                                                                                                                                                                                                                                                                                                                                                     GetAtomTypes().HA2, 0.075, 0.000, 0.000, 0.000,  0.000,  6.070,/*G*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 7.50, 3.55,  7.50,  3.70,  9.60, 2.34, 0.130, 0.150,  0.109,  0.081, -0.153,  1.150, -1.000, -0.400,  0.480,  0.000,  0.330, -0.300,  0.010,  0.080, -0.085, -0.019,  0.123, -0.048,  0.071,  0.408, -0.210,  0.113,  0.355, -0.105, -0.172, -0.100, -1.309,  0.544, -0.841, -0.122, -0.989, 185.154, 1.08,  0.774,   0.00,  29.86,    0,   0, false, 0.91))), //GLY
      HIS( AddEnum( "HISTIDINE",             /*  8 */   AATypeData( "HIS", 'H', true,  "HIS", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().CG,  GetAtomTypes().ND1, GetAtomTypes().CD2, GetAtomTypes().CE1, GetAtomTypes().NE2,  GetAtomTypes().H,    GetAtomTypes().HA,   GetAtomTypes().HB2,  GetAtomTypes().HB3,  GetAtomTypes().HD1,  GetAtomTypes().HD2,  GetAtomTypes().HE1,   GetAtomTypes().HE2),                                                                                                                                GetAtomTypes().CB,  0.025, 2.990, 0.230, 4.660,  0.130,  7.690,/*H*/    1.0,  6.5,  6.5,  6.0,  6.4,  6.00,  0.0,  6.04,  6.00,  6.6,  5.98,  6.08, 7.50, 3.55,  4.89,  6.89,  9.17, 1.82, 0.270, 0.300,  0.005,  0.081, -0.075,  0.110,  3.000, -3.200, -0.400,  0.500, -0.500,  0.100, -0.020,  0.060,  0.191, -0.110, -0.040,  0.155, -0.125,  0.001,  0.413, -0.180,  0.042,  0.061, -0.067, -0.099, -1.672,  0.538,  0.976, -0.237,  2.773, 290.040, 5.51,  9.716,  28.68,  98.73,    2,   2, true , 1.76))), //HIS
      ILE( AddEnum( "ISOLEUCINE",            /*  9 */   AATypeData( "ILE", 'I', true,  "ILE", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().CG1, GetAtomTypes().CG2, GetAtomTypes().CD1, GetAtomTypes().H,   GetAtomTypes().HA,   GetAtomTypes().HB,   GetAtomTypes().HG12, GetAtomTypes().HG13, GetAtomTypes().HG21, GetAtomTypes().HG22, GetAtomTypes().HG23, GetAtomTypes().HD11,  GetAtomTypes().HD12, GetAtomTypes().HD13),                                                                                                          GetAtomTypes().CB,  0.056, 4.190, 0.190, 4.000,  1.800,  6.040,/*I*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 7.50, 3.55,  7.48,  3.72,  9.60, 2.36, 0.300, 0.450, -0.102, -0.058,  0.204, -1.120, -3.100,  4.500,  1.380,  1.800, -1.130, -0.700, -0.280, -0.290, -0.208,  0.125,  0.157, -0.359,  0.012,  0.148,  0.144,  0.164, -0.172, -0.061,  0.170,  0.248, -3.010,  3.748,  2.775, -0.313, -1.551, 273.462, 5.11,  7.727,   0.00, 101.00,    0,   0, false, 1.93))), //ILE
      LEU( AddEnum( "LEUCINE",               /* 10 */   AATypeData( "LEU", 'L', true,  "LEU", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().CG,  GetAtomTypes().CD1, GetAtomTypes().CD2, GetAtomTypes().H,   GetAtomTypes().HA,   GetAtomTypes().HB2,  GetAtomTypes().HB3,  GetAtomTypes().HG,   GetAtomTypes().HD11, GetAtomTypes().HD12, GetAtomTypes().HD13, GetAtomTypes().HD21,  GetAtomTypes().HD22, GetAtomTypes().HD23),                                                                                                          GetAtomTypes().CB,  0.088, 2.590, 0.190, 4.000,  1.700,  6.040,/*L*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 7.50, 3.55,  7.46,  3.73,  9.60, 2.36, 0.390, 0.310, -0.164,  0.028,  0.191, -1.250, -2.800,  3.800,  1.060,  1.800, -1.180, -0.500, -0.280, -0.360, -0.187,  0.094,  0.151, -0.323, -0.089, -0.012, -0.042,  0.295,  0.024, -0.028,  0.202,  0.226, -2.199,  4.597,  3.560, -0.013, -1.913, 278.520, 4.70,  7.727,   0.00, 106.10,    0,   0, false, 1.80))), //LEU
      LYS( AddEnum( "LYSINE",                /* 11 */   AATypeData( "LYS", 'K', true,  "LYS", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().CG,  GetAtomTypes().CD,  GetAtomTypes().CE,  GetAtomTypes().NZ,  GetAtomTypes().H,    GetAtomTypes().HA,   GetAtomTypes().HB2,  GetAtomTypes().HB3,  GetAtomTypes().HG2,  GetAtomTypes().HG3,  GetAtomTypes().HD2,  GetAtomTypes().HD3,   GetAtomTypes().HE2,  GetAtomTypes().HE3,   GetAtomTypes().HZ1, GetAtomTypes().HZ2,  GetAtomTypes().HZ3),                                            GetAtomTypes().CB,  0.059, 1.890, 0.220, 4.770, -0.990,  9.990,/*K*/    1.0, 10.8, 10.0, 10.5, 10.4, 11.50, 11.2, 10.54, 10.50, 10.5, 10.00,  9.80, 7.50, 3.55, 10.00, 10.30,  8.95, 2.18, 0.320, 0.270,  0.143, -0.041, -0.078,  2.800,  8.800, -3.900, -1.500, -3.000,  1.400,  1.800,  0.320,  0.240,  0.724, -0.131, -0.220,  1.143,  0.033, -0.163,  0.402, -0.366, -0.114,  0.378, -0.085, -0.162, -0.827, -4.193, -5.659,  0.105, -0.818, 303.428, 6.57,  9.562,  27.64, 126.70,    0,   1, false, 1.82))), //LYS
      MET( AddEnum( "METHIONINE",            /* 12 */   AATypeData( "MET", 'M', true,  "MET", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().CG,  GetAtomTypes().SD,  GetAtomTypes().CE,  GetAtomTypes().H,   GetAtomTypes().HA,   GetAtomTypes().HB2,  GetAtomTypes().HB3,  GetAtomTypes().HG2,  GetAtomTypes().HG3,  GetAtomTypes().HE1,  GetAtomTypes().HE2,  GetAtomTypes().HE3),                                                                                                                                                      GetAtomTypes().CB,  0.021, 2.350, 0.220, 4.430,  1.230,  5.710,/*M*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 7.00, 3.55,  6.98,  3.68,  9.21, 2.28, 0.380, 0.320, -0.172,  0.098,  0.121, -0.670, -3.400,  1.900,  0.640,  1.300, -1.590, -0.400, -0.260, -0.190, -0.140,  0.038,  0.136, -0.305, -0.160, -0.009,  0.199,  0.471, -0.003, -0.119,  0.098,  0.170,  0.076,  3.655,  0.973, -0.219,  0.613, 291.524, 5.66,  8.892,  25.30, 102.80,    0,   0, false, 1.82))), //MET
      PHE( AddEnum( "PHENYLALANINE",         /* 13 */   AATypeData( "PHE", 'F', true,  "PHE", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().CG,  GetAtomTypes().CD1, GetAtomTypes().CD2, GetAtomTypes().CE1, GetAtomTypes().CE2,  GetAtomTypes().CZ,   GetAtomTypes().H,    GetAtomTypes().HA,   GetAtomTypes().HB2,  GetAtomTypes().HB3,  GetAtomTypes().HD1,  GetAtomTypes().HD2,   GetAtomTypes().HE1,  GetAtomTypes().HE2,   GetAtomTypes().HZ),                                                                                      GetAtomTypes().CB,  0.039, 2.940, 0.290, 5.890,  1.790,  5.670,/*F*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 7.50, 3.55,  6.96,  3.98,  9.13, 1.83, 0.300, 0.380, -0.081, -0.043,  0.149, -1.710, -3.700,  2.800,  1.190,  2.500, -2.120, -0.500, -0.410, -0.220, -0.169, -0.034,  0.294, -0.260, -0.089,  0.347,  0.014,  0.092, -0.038, -0.108, -0.055,  0.381, -2.100,  4.498,  3.397, -0.038, -0.894, 311.302, 6.16, 12.426,   0.00, 120.40,    0,   0, true , 1.77))), //PHE
      PRO( AddEnum( "PROLINE",               /* 14 */   AATypeData( "PRO", 'P', true,  "PRO", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().CG,  GetAtomTypes().CD,  GetAtomTypes().HA,  GetAtomTypes().HB2, GetAtomTypes().HB3,  GetAtomTypes().HG2,  GetAtomTypes().HG3,  GetAtomTypes().HD2,  GetAtomTypes().HD3),                                                                                                                                                                                                                     GetAtomTypes().CB,  0.046, 2.670, 0.000, 2.720,  0.720,  6.800,/*P*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 8.36, 3.55,  8.36,  3.61, 10.60, 1.99, 0.130, 0.340,  0.192,  0.550, -0.371,  0.140,  0.200, -1.600,  0.120,  0.000,  0.730,  0.300,  0.130,  0.150,  0.358, -0.103, -0.137,  0.253,  0.047,  0.115,  0.870,  0.398,  0.374, -0.005, -0.440, -0.463, -6.914, -0.757,  1.576,  0.049,  1.796, 235.409, 4.12,  5.505,   0.00,  79.84,    0,   0, false, 1.71))), //PRO
      SER( AddEnum( "SERINE",                /* 15 */   AATypeData( "SER", 'S', true,  "SER", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().OG,  GetAtomTypes().H,   GetAtomTypes().HA,  GetAtomTypes().HB2, GetAtomTypes().HB3,  GetAtomTypes().HG),                                                                                                                                                                                                                                                                                     GetAtomTypes().CB,  0.060, 1.310, 0.060, 1.600, -0.040,  5.700,/*S*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 6.93, 3.55,  6.86,  3.61,  9.15, 2.21, 0.200, 0.280,  0.089, -0.032, -0.048,  0.460, -0.600, -0.800, -0.180, -0.300,  0.520,  0.100,  0.050,  0.160,  0.050, -0.023, -0.023,  0.184,  0.125,  0.051, -0.037, -0.141,  0.043, -0.088, -0.058, -0.008, -0.995, -1.692, -3.127,  0.316,  1.191, 223.038, 2.66,  2.859,  20.23,  52.65,    1,   1, false, 1.69))), //SER
      THR( AddEnum( "THREONINE",             /* 16 */   AATypeData( "THR", 'T', true,  "THR", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().OG1, GetAtomTypes().CG2, GetAtomTypes().H,   GetAtomTypes().HA,  GetAtomTypes().HB,   GetAtomTypes().HG1,  GetAtomTypes().HG21, GetAtomTypes().HG22, GetAtomTypes().HG23),                                                                                                                                                                                                                    GetAtomTypes().CB,  0.055, 3.030, 0.110, 2.600,  0.260,  5.600,/*T*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 6.82, 3.55,  7.02,  3.57,  9.10, 2.09, 0.210, 0.360,  0.071, -0.085,  0.025,  0.250, -1.200, -0.700, -0.050,  0.400,  0.070,  0.200,  0.020, -0.080,  0.003,  0.009, -0.012,  0.055,  0.087,  0.062,  0.012, -0.169, -0.085, -0.010,  0.055,  0.042, -1.294,  0.076, -2.225,  0.411, -0.820, 243.554, 3.93,  4.694,  20.23,  72.23,    1,   1, false, 1.83))), //THR
      TRP( AddEnum( "TRYPTOPHAN",            /* 17 */   AATypeData( "TRP", 'W', true,  "TRP", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().CG,  GetAtomTypes().CD1, GetAtomTypes().CD2, GetAtomTypes().NE1, GetAtomTypes().CE2,  GetAtomTypes().CE3,  GetAtomTypes().CZ2,  GetAtomTypes().CZ3,  GetAtomTypes().CH2,  GetAtomTypes().H,    GetAtomTypes().HA,   GetAtomTypes().HB2,   GetAtomTypes().HB3,  GetAtomTypes().HD1,   GetAtomTypes().HE1, GetAtomTypes().HE3,  GetAtomTypes().HZ2, GetAtomTypes().HZ3,  GetAtomTypes().HH2),   GetAtomTypes().CB,  0.013, 3.210, 0.410, 8.080,  2.250,  5.940,/*W*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 7.50, 3.55,  7.11,  3.78,  9.39, 2.83, 0.320, 0.420, -0.028, -0.069,  0.113, -2.090, -1.900, -0.900,  0.810,  3.400, -0.510, -0.300, -0.150, -0.280,  0.025, -0.214,  0.299,  0.024, -0.232,  0.397,  0.047, -0.316,  0.115,  0.110, -0.108,  0.346, -1.870,  4.387,  2.867, -0.007,  0.444, 350.681, 7.38, 17.695,  15.79, 148.64,    0,   1, true , 1.80))), //TRP
      TYR( AddEnum( "TYROSINE",              /* 18 */   AATypeData( "TYR", 'Y', true,  "TYR", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().CG,  GetAtomTypes().CD1, GetAtomTypes().CD2, GetAtomTypes().CE1, GetAtomTypes().CE2,  GetAtomTypes().CZ,   GetAtomTypes().OH,   GetAtomTypes().H,    GetAtomTypes().HA,   GetAtomTypes().HB2,  GetAtomTypes().HB3,  GetAtomTypes().HD1,   GetAtomTypes().HD2,  GetAtomTypes().HE1,   GetAtomTypes().HE2, GetAtomTypes().HH),                                                                  GetAtomTypes().CB,  0.034, 2.940, 0.300, 6.470,  0.960,  5.660,/*Y*/   -1.0, 10.1, 10.0, 10.1, 10.0, 10.07,  0.0, 10.46, 10.00, 10.3, 10.00,  9.84, 7.50, 3.55,  9.34, 10.34,  9.11, 2.20, 0.250, 0.410,  0.112, -0.211,  0.177, -0.710,  0.700, -1.300,  0.260,  2.300, -0.210,  0.400, -0.090, -0.030, -0.038, -0.012,  0.053,  0.247, -0.015,  0.136, -0.382,  0.017, -0.114,  0.008,  0.173,  0.221, -0.445,  3.577,  1.147, -0.038, -2.688, 328.820, 6.83, 13.607,  20.23, 130.48,    1,   1, true , 1.77))), //TYR
      VAL( AddEnum( "VALINE",                /* 19 */   AATypeData( "VAL", 'V', true,  "VAL", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().CG1, GetAtomTypes().CG2, GetAtomTypes().H,   GetAtomTypes().HA,  GetAtomTypes().HB,   GetAtomTypes().HG11, GetAtomTypes().HG12, GetAtomTypes().HG13, GetAtomTypes().HG21, GetAtomTypes().HG22, GetAtomTypes().HG23),                                                                                                                                                                          GetAtomTypes().CB,  0.071, 3.670, 0.140, 3.000,  1.220,  6.020,/*V*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 7.44, 3.55,  7.44,  3.69,  9.62, 2.32, 0.270, 0.490,  0.005, -0.174,  0.243, -0.460, -2.600,  4.200,  1.080,  1.500, -1.270, -0.600, -0.170, -0.240, -0.166,  0.147,  0.066, -0.194,  0.071,  0.177, -0.044,  0.163, -0.309, -0.108,  0.311,  0.247, -2.129,  4.083,  2.297,  0.413,  5.126, 250.093, 4.27,  5.892,   0.00,  81.37,    0,   0, false, 1.87))), //VAL

      ASX( AddEnum( "ASPARAGINE_or_ASPARTIC_ACID",      AATypeData( "ASX", 'B', false, "ASP", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.052, 1.600, 0.118, 2.852, -0.698,  4.459,/*B*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,  7.46,  3.57,  9.60, 2.20, 0.233, 0.208,  0.164,  0.073, -0.181,  2.461,  7.341, -3.500, -0.849, -1.808,  0.652,  0.557,  0.289,  0.329,  0.499, -0.100, -0.185,  0.876,  0.154, -0.131,  0.156, -0.263,  0.145,  0.201, -0.205, -0.220,  0.000,  0.000,  0.000,  0.000,  0.000, 258.776, 4.16, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))), //ASX weighted average according to natural prevalence w/asn = 0.422556  asp = 0.577444
      GLX( AddEnum( "GLUTAMINE_or_GLUTAMIC_ACID",       AATypeData( "GLX", 'Z', false, "GLU", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.058, 1.560, 0.161, 3.841, -0.490,  4.006,/*Z*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,  6.96,  3.54,  9.60, 2.20, 0.399, 0.224,  0.018,  0.035, -0.046,  2.607,  6.734, -3.500, -0.779, -2.006,  0.873,  0.700,  0.189,  0.204,  0.492,  0.032, -0.276,  0.867,  0.054, -0.365,  0.100, -0.064, -0.054,  0.229,  0.128, -0.195,  0.000,  0.000,  0.000,  0.000,  0.000, 285.646, 5.22, 0.0000,  0.000,  0.000,    0,   0, false, 2.05))), //GLX weighted average according to natural prevalence w/gln = 0.35762   glu = 0.64238
      XXX( AddEnum( "ARBITRARY_AMINO_ACID",             AATypeData( "XXX", 'X', false, "XXX", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.059, 2.165, 0.147, 3.330,  0.358,  6.141,/*X*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,  7.26,  3.57,  9.60, 2.20, 0.290, 0.298,  0.024,  0.029,  0.013,  0.624,  1.524, -0.292,  0.016, -0.042,  0.023,  0.138,  0.003,  0.008,  0.110,  0.007, -0.015,  0.241,  0.032,  0.011,  0.097,  0.047,  0.026,  0.009,  0.038,  0.020,  0.000,  0.000,  0.000,  0.000,  0.000, 263.191, 4.45, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))), //XXX weighted averaged according to natural prevalence)
      UNK( AddEnum( "UNKNOWN_AMINO_ACID",               AATypeData( "UNK", 'U', false, "UNK", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.059, 2.165, 0.147, 3.330,  0.358,  6.141,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,  0.00,  0.00,  9.60, 2.20, 0.290, 0.298,  0.024,  0.029,  0.013,  0.624,  1.524, -0.292,  0.016, -0.042,  0.023,  0.138,  0.003,  0.008,  0.110,  0.007, -0.015,  0.241,  0.032,  0.011,  0.097,  0.047,  0.026,  0.009,  0.038,  0.020,  0.000,  0.000,  0.000,  0.000,  0.000, 263.191, 4.45, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))), //UNK weighted averaged according to natural prevalence)
      GAP( AddEnum( "SEQUENCE_GAP",                     AATypeData( "GAP", '-', false, "GAP"))),

      DAL( AddEnum( "D-ALANINE",                        AATypeData( "DAL", 'U', false, "ALA", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.085, 1.280, 0.050, 1.000,  0.310,  6.110,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.420, 0.230, -0.168,  0.097,  0.116,  0.500, -1.600,  1.800,  0.620,  0.500,  0.100, -0.300, -0.170, -0.150, -0.130,  0.075,  0.081, -0.207, -0.071, -0.136, -0.109,  0.215,  0.211, -0.044,  0.177,  0.153,  0.000,  0.000,  0.000,  0.000,  0.000, 209.020, 2.15, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      ALS( AddEnum( "2-AMINO-3-OXO-4-SULFO-BUTYRIC ACID",
          AATypeData( "ALS", 'U', false, "ALA", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.085, 1.280, 0.050, 1.000,  0.310,  6.110,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.420, 0.230, -0.168,  0.097,  0.116,  0.500, -1.600,  1.800,  0.620,  0.500,  0.100, -0.300, -0.170, -0.150, -0.130,  0.075,  0.081, -0.207, -0.071, -0.136, -0.109,  0.215,  0.211, -0.044,  0.177,  0.153,  0.000,  0.000,  0.000,  0.000,  0.000, 209.020, 2.15, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      ACL( AddEnum( "DEOXY-CHLOROMETHYL-ARGININE",      AATypeData( "ACL", 'U', false, "ARG", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.051, 2.340, 0.290, 6.130, -1.010, 10.740,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.360, 0.250,  0.076, -0.033, -0.036,  1.810, 12.300, -4.500, -2.530, -3.000,  1.910,  1.400,  0.370,  0.320,  0.400, -0.110, -0.150,  0.978, -0.093, -0.149,  0.056, -0.084, -0.069, -0.002, -0.018, -0.068,  0.000,  0.000,  0.000,  0.000,  0.000, 335.732, 7.36, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      CME( AddEnum( "S,S-(2-HYDROXYETHYL)THIOCYSTEINE", AATypeData( "CME", 'U', false, "CYS", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.013, 1.770, 0.130, 2.430,  1.540,  6.350,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.170, 0.410, -0.173,  0.482, -0.114, -0.020, -2.000,  2.500,  0.290,  1.000, -1.420, -0.900, -0.060, -0.150, -0.070,  0.560, -0.220, -0.347,  0.382, -0.375,  1.925,  1.454, -0.058, -0.169,  0.463, -0.325,  0.000,  0.000,  0.000,  0.000,  0.000, 240.500, 3.71, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      CSE( AddEnum( "SELENOCYSTEINE",                   AATypeData( "CSE", 'U', false, "CYS", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.013, 1.770, 0.130, 2.430,  1.540,  6.350,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.170, 0.410, -0.030,  0.245, -0.148, -0.020, -2.000,  2.500,  0.290,  1.000, -1.420, -0.900, -0.060, -0.150,  0.040,  0.369, -0.250, -0.125,  0.240, -0.285,  1.381,  1.205, -0.159, -0.228,  0.281, -0.353,  0.000,  0.000,  0.000,  0.000,  0.000, 240.500, 3.71, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      CSD( AddEnum( "3-SULFINOALANINE",                 AATypeData( "CSD", 'U', false, "CYS", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.013, 1.770, 0.130, 2.430,  1.540,  6.350,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.170, 0.410, -0.173,  0.482, -0.114, -0.020, -2.000,  2.500,  0.290,  1.000, -1.420, -0.900, -0.060, -0.150, -0.070,  0.560, -0.220, -0.347,  0.382, -0.375,  1.925,  1.454, -0.058, -0.169,  0.463, -0.325,  0.000,  0.000,  0.000,  0.000,  0.000, 240.500, 3.71, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      CSO( AddEnum( "S-HYDROXYCYSTEINE",                AATypeData( "CSO", 'U', false, "CYS", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.013, 1.770, 0.130, 2.430,  1.540,  6.350,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.170, 0.410, -0.173,  0.482, -0.114, -0.020, -2.000,  2.500,  0.290,  1.000, -1.420, -0.900, -0.060, -0.150, -0.070,  0.560, -0.220, -0.347,  0.382, -0.375,  1.925,  1.454, -0.058, -0.169,  0.463, -0.325,  0.000,  0.000,  0.000,  0.000,  0.000, 240.500, 3.71, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      CSW( AddEnum( "CYSTEINE-S-DIOXIDE",               AATypeData( "CSW", 'U', false, "CYS", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.013, 1.770, 0.130, 2.430,  1.540,  6.350,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.170, 0.410, -0.173,  0.482, -0.114, -0.020, -2.000,  2.500,  0.290,  1.000, -1.420, -0.900, -0.060, -0.150, -0.070,  0.560, -0.220, -0.347,  0.382, -0.375,  1.925,  1.454, -0.058, -0.169,  0.463, -0.325,  0.000,  0.000,  0.000,  0.000,  0.000, 240.500, 3.71, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      CYG( AddEnum( "2-AMINO-4-(AMINO-3-OXO-PROPYLSULFANYLCARBONYL)-BUTYRIC ACID",
          AATypeData( "CYG", 'U', false, "CYS", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.013, 1.770, 0.130, 2.430,  1.540,  6.350,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.170, 0.410, -0.173,  0.482, -0.114, -0.020, -2.000,  2.500,  0.290,  1.000, -1.420, -0.900, -0.060, -0.150, -0.070,  0.560, -0.220, -0.347,  0.382, -0.375,  1.925,  1.454, -0.058, -0.169,  0.463, -0.325,  0.000,  0.000,  0.000,  0.000,  0.000, 240.500, 3.71, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      OCS( AddEnum( "CYSTEINESULFONIC_ACID",            AATypeData( "OCS", 'U', false, "CYS", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.013, 1.770, 0.130, 2.430,  1.540,  6.350,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.170, 0.410, -0.173,  0.482, -0.114, -0.020, -2.000,  2.500,  0.290,  1.000, -1.420, -0.900, -0.060, -0.150, -0.070,  0.560, -0.220, -0.347,  0.382, -0.375,  1.925,  1.454, -0.058, -0.169,  0.463, -0.325,  0.000,  0.000,  0.000,  0.000,  0.000, 240.500, 3.71, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      SC2( AddEnum( "N-ACETYL-L-CYSTEINE",              AATypeData( "SC2", 'U', false, "CYS", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.013, 1.770, 0.130, 2.430,  1.540,  6.350,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.170, 0.410, -0.173,  0.482, -0.114, -0.020, -2.000,  2.500,  0.290,  1.000, -1.420, -0.900, -0.060, -0.150, -0.070,  0.560, -0.220, -0.347,  0.382, -0.375,  1.925,  1.454, -0.058, -0.169,  0.463, -0.325,  0.000,  0.000,  0.000,  0.000,  0.000, 240.500, 3.71, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      CGU( AddEnum( "GAMMA-CARBOXY-GLUTAMIC_ACID",      AATypeData( "CGU", 'U', false, "GLU", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.069, 1.560, 0.150, 3.780, -0.640,  3.090,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.420, 0.210,  0.029,  0.058, -0.078,  3.630,  8.200, -3.500, -0.740, -3.000,  0.830,  0.700,  0.150,  0.300,  0.470,  0.100, -0.310,  0.828,  0.192, -0.352,  0.235, -0.025, -0.034, -0.184,  0.196, -0.175,  0.000,  0.000,  0.000,  0.000,  0.000, 285.025, 5.05, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      PCA( AddEnum( "PYROGLUTAMIC_ACID",                AATypeData( "PCA", 'U', false, "GLU", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.069, 1.560, 0.150, 3.780, -0.640,  3.090,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.420, 0.210,  0.029,  0.058, -0.078,  3.630,  8.200, -3.500, -0.740, -3.000,  0.830,  0.700,  0.150,  0.300,  0.470,  0.100, -0.310,  0.828,  0.192, -0.352,  0.235, -0.025, -0.034, -0.184,  0.196, -0.175,  0.000,  0.000,  0.000,  0.000,  0.000, 285.025, 5.05, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      KCX( AddEnum( "LYSINE_NZ-CARBOXYLIC ACID",        AATypeData( "KCX", 'U', false, "LYS", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.059, 1.890, 0.220, 4.770, -0.990,  9.990,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.320, 0.270,  0.152, -0.043, -0.083,  2.800,  8.800, -3.900, -1.500, -3.000,  1.400,  1.800,  0.320,  0.240,  0.690, -0.130, -0.220,  1.253,  0.102, -0.129,  0.425, -0.360, -0.100, -0.025,  0.025, -0.137,  0.000,  0.000,  0.000,  0.000,  0.000, 303.428, 6.57, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      LLP( AddEnum( "2-LYSINE(3-HYDROXY-2-METHYL-5-PHOSPHONOOXYMETHYL-PYRIDIN-4-YLMETHANE)",
          AATypeData( "LLP", 'U', false, "LYS", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.059, 1.890, 0.220, 4.770, -0.990,  9.990,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.320, 0.270,  0.152, -0.043, -0.083,  2.800,  8.800, -3.900, -1.500, -3.000,  1.400,  1.800,  0.320,  0.240,  0.690, -0.130, -0.220,  1.253,  0.102, -0.129,  0.425, -0.360, -0.100, -0.025,  0.025, -0.137,  0.000,  0.000,  0.000,  0.000,  0.000, 303.428, 6.57, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      M3L( AddEnum( "N-TRIMETHYLLYSINE",                AATypeData( "M3L", 'U', false, "LYS", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.059, 1.890, 0.220, 4.770, -0.990,  9.990,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.320, 0.270,  0.143, -0.041, -0.078,  2.800,  8.800, -3.900, -1.500, -3.000,  1.400,  1.800,  0.320,  0.240,  0.724, -0.131, -0.220,  1.143,  0.033, -0.163,  0.402, -0.366, -0.114,  0.378, -0.085, -0.162,  0.000,  0.000,  0.000,  0.000,  0.000, 303.428, 6.57, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      MLY( AddEnum( "N-DIMETHYL-LYSINE",                AATypeData( "MLY", 'U', false, "LYS", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.059, 1.890, 0.220, 4.770, -0.990,  9.990,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.320, 0.270,  0.152, -0.043, -0.083,  2.800,  8.800, -3.900, -1.500, -3.000,  1.400,  1.800,  0.320,  0.240,  0.690, -0.130, -0.220,  1.253,  0.102, -0.129,  0.425, -0.360, -0.100, -0.025,  0.025, -0.137,  0.000,  0.000,  0.000,  0.000,  0.000, 303.428, 6.57, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      LYR( AddEnum( "N~6~-[(2Z,4E,6E,8E)-3,7-DIMETHYL-9-(2,6,6-TRIMETHYLCYCLOHEX-1-EN-1-YL)NONA-2,4,6,8-TETRAENYL]LYSINE",
          AATypeData( "LYR", 'U', false, "LYS", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.059, 1.890, 0.220, 4.770, -0.990,  9.990,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.320, 0.270,  0.152, -0.043, -0.083,  2.800,  8.800, -3.900, -1.500, -3.000,  1.400,  1.800,  0.320,  0.240,  0.690, -0.130, -0.220,  1.253,  0.102, -0.129,  0.425, -0.360, -0.100, -0.025,  0.025, -0.137,  0.000,  0.000,  0.000,  0.000,  0.000, 303.428, 6.57, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      CXM( AddEnum( "N-CARBOXYMETHIONINE",              AATypeData( "CXM", 'U', false, "MET", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.021, 2.350, 0.220, 4.430,  1.230,  5.710,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.380, 0.320, -0.174,  0.074,  0.151, -0.670, -3.400,  1.900,  0.640,  1.300, -1.590, -0.400, -0.260, -0.190, -0.120,  0.010,  0.130, -0.255, -0.183,  0.037,  0.139,  0.441, -0.045, -0.221,  0.131,  0.354,  0.000,  0.000,  0.000,  0.000,  0.000, 291.524, 5.66, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      FME( AddEnum( "N-FORMYLMETHIONINE",               AATypeData( "FME", 'U', false, "MET", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.021, 2.350, 0.220, 4.430,  1.230,  5.710,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.380, 0.320, -0.174,  0.074,  0.151, -0.670, -3.400,  1.900,  0.640,  1.300, -1.590, -0.400, -0.260, -0.190, -0.120,  0.010,  0.130, -0.255, -0.183,  0.037,  0.139,  0.441, -0.045, -0.221,  0.131,  0.354,  0.000,  0.000,  0.000,  0.000,  0.000, 291.524, 5.66, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      MSE( AddEnum( "SELENO_METHIONINE",                AATypeData( "MSE", 'U', false, "MET", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().CG,  GetAtomTypes().SE,  GetAtomTypes().CE,  GetAtomTypes().H,   GetAtomTypes().HA,   GetAtomTypes().HB2,  GetAtomTypes().HB3,  GetAtomTypes().HG2,  GetAtomTypes().HG3,  GetAtomTypes().HE1,  GetAtomTypes().HE2,  GetAtomTypes().HE3),                                                                                                                                                      GetAtomTypes().CB,  0.021, 2.350, 0.220, 4.430,  1.230,  5.710,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.380, 0.320, -0.174,  0.074,  0.151, -0.670, -3.400,  1.900,  0.640,  1.300, -1.590, -0.400, -0.260, -0.190, -0.120,  0.010,  0.130, -0.255, -0.183,  0.037,  0.139,  0.441, -0.045, -0.221,  0.131,  0.354,  0.000,  0.000,  0.000,  0.000,  0.000, 291.524, 5.66, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      DPN( AddEnum( "D-PHENYLALANINE",                  AATypeData( "DPN", 'U', false, "PHE", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.039, 2.940, 0.290, 5.890,  1.790,  5.670,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.300, 0.380, -0.104, -0.016,  0.148, -1.710, -3.700,  2.800,  1.190,  2.500, -2.120, -0.500, -0.410, -0.220, -0.190, -0.020,  0.340, -0.342, -0.101,  0.371,  0.031, -0.015, -0.011,  0.097, -0.084,  0.416,  0.000,  0.000,  0.000,  0.000,  0.000, 311.302, 6.16, 0.0000,  0.000,  0.000,    0,   0, true , 2.00))),
      S1H( AddEnum( "1-HEXADECANOSULFONYL-O-L-SERINE",  AATypeData( "S1H", 'U', false, "SER", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.060, 1.310, 0.060, 1.600, -0.040,  5.700,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.200, 0.280,  0.097, -0.018, -0.068,  0.460, -0.600, -0.800, -0.180, -0.300,  0.520,  0.100,  0.050,  0.160,  0.020,  0.020, -0.040,  0.160,  0.216,  0.000, -0.032, -0.076,  0.038, -0.094, -0.097, -0.033,  0.000,  0.000,  0.000,  0.000,  0.000, 223.038, 2.66, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      SAC( AddEnum( "N-ACETYL-SERINE",                  AATypeData( "SAC", 'U', false, "SER", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.060, 1.310, 0.060, 1.600, -0.040,  5.700,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.200, 0.280,  0.097, -0.018, -0.068,  0.460, -0.600, -0.800, -0.180, -0.300,  0.520,  0.100,  0.050,  0.160,  0.020,  0.020, -0.040,  0.160,  0.216,  0.000, -0.032, -0.076,  0.038, -0.094, -0.097, -0.033,  0.000,  0.000,  0.000,  0.000,  0.000, 223.038, 2.66, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      SEP( AddEnum( "PHOSPHOSERINE",                    AATypeData( "SEP", 'U', false, "SER", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.060, 1.310, 0.060, 1.600, -0.040,  5.700,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.200, 0.280,  0.097, -0.018, -0.068,  0.460, -0.600, -0.800, -0.180, -0.300,  0.520,  0.100,  0.050,  0.160,  0.020,  0.020, -0.040,  0.160,  0.216,  0.000, -0.032, -0.076,  0.038, -0.094, -0.097, -0.033,  0.000,  0.000,  0.000,  0.000,  0.000, 223.038, 2.66, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      TPO( AddEnum( "PHOSPHOTHREONINE",                 AATypeData( "TPO", 'U', false, "THR", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.055, 3.030, 0.110, 2.600,  0.260,  5.600,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.210, 0.360,  0.104, -0.100,  0.015,  0.250, -1.200, -0.700, -0.050,  0.400,  0.070,  0.200,  0.020, -0.080, -0.010,  0.020,  0.000,  0.061,  0.136,  0.128, -0.020, -0.230, -0.092,  0.010,  0.098,  0.008,  0.000,  0.000,  0.000,  0.000,  0.000, 243.554, 3.93, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      TYS( AddEnum( "O-SULFO-L-TYROSINE",               AATypeData( "TYS", 'U', false, "TYR", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.034, 2.940, 0.300, 6.470,  0.960,  5.660,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.250, 0.410,  0.140, -0.242,  0.210, -0.710,  0.700, -1.300,  0.260,  2.300, -0.210,  0.400, -0.090, -0.030, -0.070,  0.040,  0.030,  0.259, -0.010,  0.093, -0.423,  0.005, -0.182,  0.358,  0.176,  0.141,  0.000,  0.000,  0.000,  0.000,  0.000, 328.820, 6.83, 0.0000,  0.000,  0.000,    0,   0, true , 2.00))),
      R1A( AddEnum( "methanesulfonothioate",            AATypeData( "R1A", 'U', false, "CYS", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB, GetAtomTypes().SG, GetAtomTypes().SD, GetAtomTypes().CE, GetAtomTypes().C2, GetAtomTypes().C3, GetAtomTypes().C4, GetAtomTypes().C5, GetAtomTypes().C6, GetAtomTypes().C7, GetAtomTypes().C8, GetAtomTypes().C9, GetAtomTypes().N1, GetAtomTypes().O1),                                                                                                                                                      GetAtomTypes().CB,  0.000, 0.000, 0.000, 0.000,  0.000,  0.000,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.000, 0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, 000.000, 0.00, 0.0000,  0.000,  0.000,    0,   0, false, 2.00))),
      IAS( AddEnum( "BETA-L-ASPARTIC_ACID",             AATypeData( "IAS", 'U', false, "ASP", storage::Set< AtomType>::Create( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().C, GetAtomTypes().O, GetAtomTypes().CB),                                                                                                                                                                                                                                                                                                                                                                                                             GetAtomTypes().CB,  0.058, 1.600, 0.110, 2.780, -0.770,  2.950,/*U*/    0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.00,  0.00,  0.0,  0.00,  0.00, 0.00, 0.00,   0.0,   0.0,  9.60, 2.20, 0.250, 0.200,  0.169,  0.067, -0.180,  3.640,  9.200, -3.500, -0.900, -3.000,  0.780,  0.600,  0.370,  0.410,  0.574, -0.068, -0.237,  0.963,  0.227, -0.193,  0.237, -0.295,  0.074,  0.266, -0.176, -0.253,  0.000,  0.000,  0.000,  0.000,  0.000, 257.993, 4.00, 0.0000,  0.000,  0.000,    0,   0, false, 2.00)))
    {
      // @see http://www.bmrb.wisc.edu/ref_info/atom_nom.tbl
      // add mapping of atom types to pdb atom names
      ALA->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB1 , "1HB ");
      ALA->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB2 , "2HB ");
      ALA->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB3 , "3HB ");

      ARG->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB2 , "1HB ");
      ARG->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB3 , "2HB ");
      ARG->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG2 , "1HG ");
      ARG->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG3 , "2HG ");
      ARG->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HD2 , "1HD ");
      ARG->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HD3 , "2HD ");
      ARG->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HH11, "1HH1");
      ARG->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HH12, "2HH1");
      ARG->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HH21, "1HH2");
      ARG->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HH22, "2HH2");

      ASP->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB2 , "1HB ");
      ASP->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB3 , "2HB ");

      ASN->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB2 , "1HB ");
      ASN->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB3 , "2HB ");
      ASN->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HD21, "1HD2");
      ASN->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HD22, "2HD2");

      CYS->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB2 , "1HB ");
      CYS->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB3 , "2HB ");

      GLU->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB2 , "1HB ");
      GLU->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB3 , "2HB ");
      GLU->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG2 , "1HG ");
      GLU->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG3 , "2HG ");

      GLN->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB2 , "1HB ");
      GLN->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB3 , "2HB ");
      GLN->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG2 , "1HG ");
      GLN->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG3 , "2HG ");
      GLN->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HE21, "1HE2");
      GLN->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HE22, "2HE2");

      GLY->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HA2 , "1HA ");
      GLY->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HA3 , "2HA ");

      HIS->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB2 , "1HB ");
      HIS->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB3 , "2HB ");

      ILE->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG12, "1HG1");
      ILE->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG13, "2HG1");
      ILE->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG21, "1HG2");
      ILE->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG22, "2HG2");
      ILE->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG23, "3HG2");
      ILE->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HD11, "1HD1");
      ILE->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HD12, "2HD1");
      ILE->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HD13, "3HD1");

      LEU->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB2 , "1HB ");
      LEU->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB3 , "2HB ");
      LEU->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HD11, "1HD1");
      LEU->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HD12, "2HD1");
      LEU->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HD13, "3HD1");
      LEU->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HD21, "1HD2");
      LEU->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HD22, "2HD2");
      LEU->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HD23, "3HD2");

      LYS->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB2 , "1HB ");
      LYS->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB3 , "2HB ");
      LYS->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG2 , "1HG ");
      LYS->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG3 , "2HG ");
      LYS->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HD2 , "1HD ");
      LYS->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HD3 , "2HD ");
      LYS->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HE2 , "1HE ");
      LYS->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HE3 , "2HE ");
      LYS->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HZ1 , "1HZ ");
      LYS->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HZ2 , "2HZ ");
      LYS->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HZ3 , "3HZ ");

      MET->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB2 , "1HB ");
      MET->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB3 , "2HB ");
      MET->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG2 , "1HG ");
      MET->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG3 , "2HG ");
      MET->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HE1 , "1HE ");
      MET->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HE2 , "2HE ");
      MET->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HE3 , "3HE ");

      PHE->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB2 , "1HB ");
      PHE->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB3 , "2HB ");

      PRO->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB2 , "1HB ");
      PRO->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB3 , "2HB ");
      PRO->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG2 , "1HG ");
      PRO->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG3 , "2HG ");
      PRO->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HD2 , "1HD ");
      PRO->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HD3 , "2HD ");

      SER->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB2 , "1HB ");
      SER->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB3 , "2HB ");

      THR->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG21, "1HG2");
      THR->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG22, "2HG2");
      THR->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG23, "3HG2");

      TRP->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB2 , "1HB ");
      TRP->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB3 , "2HB ");

      TYR->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB2 , "1HB ");
      TYR->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HB3 , "2HB ");

      VAL->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG11, "1HG1");
      VAL->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG12, "2HG1");
      VAL->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG13, "3HG1");
      VAL->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG21, "1HG2");
      VAL->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG22, "2HG2");
      VAL->AddAtomTypePDBAtomNameMapping( GetAtomTypes().HG23, "3HG2");

      // set pka for commonly charged hydrogens
      ARG->SetSideChainIonizableHType( GetAtomTypes().HD22);
      HIS->SetSideChainIonizableHType( GetAtomTypes().HD1);
      LYS->SetSideChainIonizableHType( GetAtomTypes().HZ3);
      ASP->SetSideChainIonizableHType( GetAtomTypes().HD2);
      GLU->SetSideChainIonizableHType( GetAtomTypes().HE2);
      CYS->SetSideChainIonizableHType( GetAtomTypes().HG);
      TYR->SetSideChainIonizableHType( GetAtomTypes().HH);

      // set the dihedral angle atoms
      MET->SetDihedralAngles
      (
        storage::Map< ChiAngle::ChiEnum, storage::VectorND< 4, AtomType> >::Create
        (
          std::pair< ChiAngle::ChiEnum, storage::VectorND< 4, AtomType> >
          (
            ChiAngle::e_One,
            storage::VectorND< 4, AtomType>( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().CB, GetAtomTypes().CG)
          ),
          std::pair< ChiAngle::ChiEnum, storage::VectorND< 4, AtomType> >
          (
            ChiAngle::e_Two,
            storage::VectorND< 4, AtomType>( GetAtomTypes().CA, GetAtomTypes().CB, GetAtomTypes().CG, GetAtomTypes().SD)
          ),
          std::pair< ChiAngle::ChiEnum, storage::VectorND< 4, AtomType> >
          (
            ChiAngle::e_Three,
            storage::VectorND< 4, AtomType>( GetAtomTypes().CB, GetAtomTypes().CG, GetAtomTypes().SD, GetAtomTypes().CE)
          )
        )
      );

      R1A->SetDihedralAngles
      (
        storage::Map< ChiAngle::ChiEnum, storage::VectorND< 4, AtomType> >::Create
        (
          std::make_pair
          (
            ChiAngle::e_One,
            storage::VectorND< 4, AtomType>( GetAtomTypes().N, GetAtomTypes().CA, GetAtomTypes().CB, GetAtomTypes().SG)
          ),
          std::pair< ChiAngle::ChiEnum, storage::VectorND< 4, AtomType> >
          (
            ChiAngle::e_Two,
            storage::VectorND< 4, AtomType>( GetAtomTypes().CA, GetAtomTypes().CB, GetAtomTypes().SG, GetAtomTypes().SD)
          ),
          std::pair< ChiAngle::ChiEnum, storage::VectorND< 4, AtomType> >
          (
            ChiAngle::e_Three,
            storage::VectorND< 4, AtomType>( GetAtomTypes().CB, GetAtomTypes().SG, GetAtomTypes().SD, GetAtomTypes().CE)
          ),
          std::pair< ChiAngle::ChiEnum, storage::VectorND< 4, AtomType> >
          (
            ChiAngle::e_Four,
            storage::VectorND< 4, AtomType>( GetAtomTypes().SG, GetAtomTypes().SD, GetAtomTypes().CE, GetAtomTypes().C3)
          ),
          std::pair< ChiAngle::ChiEnum, storage::VectorND< 4, AtomType> >
          (
            ChiAngle::e_Five,
            storage::VectorND< 4, AtomType>( GetAtomTypes().SD, GetAtomTypes().CE, GetAtomTypes().C3, GetAtomTypes().C4)
          )
        )
      );

      // determine the parent types
      for( iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
      {
        ( *itr)->DetermineParentType( *this);
      }
      e_Undefined->DetermineParentType( *this);
    }

    //! @brief function to deduce AAType from three letter code of an amino acid
    //! @param THREE_LETTER_CODE three letter code descriptor for amino acid of interest
    //! @return AAType specified by given THREE_LETTER_CODE
    const AAType &AATypes::AATypeFromThreeLetterCode( const std::string &THREE_LETTER_CODE) const
    {
      // assert the provided string has correct size
      if( THREE_LETTER_CODE.size() != 3)
      {
        BCL_MessageVrb( "three letter code should consist of three letters: " + THREE_LETTER_CODE);
        return e_Undefined;
      }

      // iterate over all AATypes
      for
      (
        const_iterator type_itr( GetAATypes().Begin()), type_itr_end( AATypes::End());
        type_itr != type_itr_end;
        ++type_itr
      )
      {
        // if THREE_LETTER_CODE matches this AATypes's three letter code
        if( ( *type_itr)->GetThreeLetterCode() == THREE_LETTER_CODE)
        {
          // return this AAType
          return *type_itr;
        }
      }

      // if no match was found return undefined AAType
      return GetAATypes().e_Undefined;
    }

    //! @brief return if two amino acids have same parent
    //! this function returns true if two amino acids denoted by there three letter code have the same parent.
    //! MSE == MET would return true
    //! @param THREE_LETTER_CODE_LHS three letter code descriptor for amino acid of interest lhs
    //! @param THREE_LETTER_CODE_RHS three letter code descriptor for amino acid of interest rhs
    //! @return true if the two amino acids have same parent
    bool AATypes::HaveSameParent( const std::string &THREE_LETTER_CODE_LHS, const std::string &THREE_LETTER_CODE_RHS) const
    {
      if( THREE_LETTER_CODE_LHS == THREE_LETTER_CODE_RHS)
      {
        return true;
      }

      const AAType &lhs( AATypeFromThreeLetterCode( THREE_LETTER_CODE_LHS));
      const AAType &rhs( AATypeFromThreeLetterCode( THREE_LETTER_CODE_RHS));
      if( !lhs.IsDefined() || !rhs.IsDefined())
      {
        return false;
      }
      return lhs->GetParentType() == rhs->GetParentType();
    }

    namespace
    {
      // anonymous namespace to prevent unwanted symbol export
      //! @brief create a vector of 256 values that return, for a given character, the corresponding aa
      //! @return a vector of 256 values that return, for a given character, the corresponding aa
      storage::Vector< AAType> CreateAATypeFromOneLetterCodeVector()
      {
        storage::Vector< AAType> aatypes( 256, GetAATypes().e_Undefined);
        // iterate over all AATypes
        for
        (
          AATypes::const_iterator type_itr( GetAATypes().Begin()), type_itr_end( GetAATypes().End());
          type_itr != type_itr_end;
          ++type_itr
        )
        {
          // check that the aatypes vector was not already set; typically this happens for unnatural types,
          // usually denoted with U as the one letter code, but we want to match
          if( !aatypes( size_t( ( *type_itr)->GetOneLetterCode())).IsDefined())
          {
            aatypes( size_t( ( *type_itr)->GetOneLetterCode())) = *type_itr;
          }
        }
        return aatypes;
      }
    }

    //! @brief function to deduce AAType from one letter code of an amino acid
    //! @param ONE_LETTER_CODE one letter code descriptor for amino acid of interest
    //! @return AAType specified by the given ONE_LETTER_CODE
    const AAType &AATypes::AATypeFromOneLetterCode( const char ONE_LETTER_CODE) const
    {
      static const storage::Vector< AAType> s_aa_type_one_letter_code( CreateAATypeFromOneLetterCodeVector());
      // return the type from the vector
      return size_t( ONE_LETTER_CODE) < s_aa_type_one_letter_code.GetSize()
             ? s_aa_type_one_letter_code( size_t( ONE_LETTER_CODE))
             : GetAATypes().e_Undefined;
    }

    //! @brief gives the 20 natural amino acid types
    //! @return set of the 20 natural amino acid types
    const storage::Set< AAType> &AATypes::GetNaturalAATypes() const
    {
      static storage::Set< AAType> s_natural_aa_types;

      if( s_natural_aa_types.IsEmpty())
      {
        // iterate over AATypes
        for
        (
          const_iterator type_itr( Begin()), type_itr_end( ASX.GetIterator());
          type_itr != type_itr_end;
          ++type_itr
        )
        {
          s_natural_aa_types.Insert( *type_itr);
        }
        BCL_Assert
        (
          20 == s_natural_aa_types.GetSize(),
          "number natural residue types should be " + util::Format()( 20)
          + " but is " + util::Format()( s_natural_aa_types.GetSize())
        );
      }

      return s_natural_aa_types;
    }

    //! @brief construct on access function for all AATypes
    //! @return reference to only instance of AATypes enum
    const AATypes &GetAATypes()
    {
      return AATypes::GetEnums();
    }

  } // namespace biol

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< biol::AATypeData, biol::AATypes>;

  } // namespace util
} // namespace bcl

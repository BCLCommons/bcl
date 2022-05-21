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
#include "score/bcl_score_aa_assignment_blosum.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_base.h"
#include "linal/bcl_linal_vector_const_reference.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    //! @brief TableType as string
    //! @param TABLE_TYPE the TableType
    //! @return the string for the TableType
    const std::string &AAAssignmentBLOSUM::GetTableDescriptor( const TableType &TABLE_TYPE)
    {
      static const std::string s_descriptors[] =
      {
        "BLOSUM_90",
        "BLOSUM_80",
        "BLOSUM_62",
        "BLOSUM_45",
        GetStaticClassName< TableType>()
      };

      return s_descriptors[ TABLE_TYPE];
    }

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AAAssignmentBLOSUM::s_Instance
    (
      GetObjectInstances().AddInstance( new AAAssignmentBLOSUM())
    );

    //! Pam and Blosum scoring matrices for amino acid replacement. Given is
    //! the log of probability of replacing aa i with aa j divided by freqency of aa i
    //! NOT mutiplied by 10 (as usually done). Note that with increasing sequence
    //! diversity move from Pam100 (45% sequence identity) to Pam250 (20% sequence identity)
    //! or from Blosum90 to Blosum45.
    //!
    //! Differences between PAM and BLOSSUM:
    //! PAM matrices are based on an explicit evolutionary model (that is, replacements
    //! are counted on the branches of a phylogenetic tree), whereas the Blosum matrices
    //! are based on an implicit rather than explicit model of evolution.
    //! The sequence variability in the alignments used to count replacements. The PAM
    //! matrices are based on mutations observed throughout a global alignment, this
    //! includes both highly conserved and highly mutable regions. The Blosum matrices
    //! are based only on highly conserved regions in series of alignments forbidden to
    //! contain gaps.
    //! The method used to count the replacements is different, unlike the PAM matrix,
    //! the Blosum procedure uses groups of sequences within which not all mutations are
    //! counted the same.
    //!
    //! Equivalent PAM and Blossum matrices: The following matrices are roughly equivalent...
    //! PAM100 <==> BLOSUM90
    //! PAM120 <==> BLOSUM80
    //! PAM160 <==> BLOSUM62
    //! PAM250 <==> BLOSUM45
    //!
    //! Generally speaking...
    //! The BLOSUM matrices are best for detecting local alignments.
    //! The BLOSUM62 matrix is the best for detecting the majority of weak protein similarities.
    //! The BLOSUM45 matrix is the best for detecting long and weak alignments.
    //!
    //! BLOSUM (Blocks Substitution Matrix):
    //! The BLOSUM matrices, also used for protein database search scoring (the default in
    //! blastp), are divided into statistical significance degrees which, in a way, are
    //! reminiscent of PAM distances. BLOSSUM matrices are most sensitive for local
    //! alignment of related sequences. The BLOSUM matrices are therefore ideal when trying
    //! to identify an unknown nucleotide sequence.
    const double AAAssignmentBLOSUM::s_BLOSUMTable[][ biol::AATypes::s_NumberStandardAATypes + 4][ biol::AATypes::s_NumberStandardAATypes + 4] =
    {
      //BLOSUM_90
      {
//        ALA  ARG  ASN  ASP  CYS  GLN  GLU  GLY  HIS  ILE  LEU  LYS  MET  PHE  PRO  SER  THR  TRP  TYR  VAL  ASX  GLX  XXX  GAP
        { 0.5,-0.2,-0.2,-0.3,-0.1,-0.1,-0.1, 0.0,-0.2,-0.2,-0.2,-0.1,-0.2,-0.3,-0.1, 0.1, 0.0,-0.4,-0.3,-0.1,-0.2,-0.1,-0.1,-0.6}, //ALA
        {-0.2, 0.6,-0.1,-0.3,-0.5, 0.1,-0.1,-0.3, 0.0,-0.4,-0.3, 0.2,-0.2,-0.4,-0.3,-0.1,-0.2,-0.4,-0.3,-0.3,-0.2, 0.0,-0.2,-0.6}, //ARG
        {-0.2,-0.1, 0.7, 0.1,-0.4, 0.0,-0.1,-0.1, 0.0,-0.4,-0.4, 0.0,-0.3,-0.4,-0.3, 0.0, 0.0,-0.5,-0.3,-0.4, 0.4,-0.1,-0.2,-0.6}, //ASN
        {-0.3,-0.3, 0.1, 0.7,-0.5,-0.1, 0.1,-0.2,-0.2,-0.5,-0.5,-0.1,-0.4,-0.5,-0.3,-0.1,-0.2,-0.6,-0.4,-0.5, 0.4, 0.0,-0.2,-0.6}, //ASP
        {-0.1,-0.5,-0.4,-0.5, 0.9,-0.4,-0.6,-0.4,-0.5,-0.2,-0.2,-0.4,-0.2,-0.3,-0.4,-0.2,-0.2,-0.4,-0.4,-0.2,-0.4,-0.5,-0.3,-0.6}, //CYS
        {-0.1, 0.1, 0.0,-0.1,-0.4, 0.7, 0.2,-0.3, 0.1,-0.4,-0.3, 0.1, 0.0,-0.4,-0.2,-0.1,-0.1,-0.3,-0.3,-0.3,-0.1, 0.4,-0.1,-0.6}, //GLN
        {-0.1,-0.1,-0.1, 0.1,-0.6, 0.2, 0.6,-0.3,-0.1,-0.4,-0.4, 0.0,-0.3,-0.5,-0.2,-0.1,-0.1,-0.5,-0.4,-0.3, 0.0, 0.4,-0.2,-0.6}, //GLU
        { 0.0,-0.3,-0.1,-0.2,-0.4,-0.3,-0.3, 0.6,-0.3,-0.5,-0.5,-0.2,-0.4,-0.5,-0.3,-0.1,-0.3,-0.4,-0.5,-0.5,-0.2,-0.3,-0.2,-0.6}, //GLY
        {-0.2, 0.0, 0.0,-0.2,-0.5, 0.1,-0.1,-0.3, 0.8,-0.4,-0.4,-0.1,-0.3,-0.2,-0.3,-0.2,-0.2,-0.3, 0.1,-0.4,-0.1, 0.0,-0.2,-0.6}, //HIS
        {-0.2,-0.4,-0.4,-0.5,-0.2,-0.4,-0.4,-0.5,-0.4, 0.5, 0.1,-0.4, 0.1,-0.1,-0.4,-0.3,-0.1,-0.4,-0.2, 0.3,-0.5,-0.4,-0.2,-0.6}, //ILE
        {-0.2,-0.3,-0.4,-0.5,-0.2,-0.3,-0.4,-0.5,-0.4, 0.1, 0.5,-0.3, 0.2, 0.0,-0.4,-0.3,-0.2,-0.3,-0.2, 0.0,-0.5,-0.4,-0.2,-0.6}, //LEU
        {-0.1, 0.2, 0.0,-0.1,-0.4, 0.1, 0.0,-0.2,-0.1,-0.4,-0.3, 0.6,-0.2,-0.4,-0.2,-0.1,-0.1,-0.5,-0.3,-0.3,-0.1, 0.1,-0.1,-0.6}, //LYS
        {-0.2,-0.2,-0.3,-0.4,-0.2, 0.0,-0.3,-0.4,-0.3, 0.1, 0.2,-0.2, 0.7,-0.1,-0.3,-0.2,-0.1,-0.2,-0.2, 0.0,-0.4,-0.2,-0.1,-0.6}, //MET
        {-0.3,-0.4,-0.4,-0.5,-0.3,-0.4,-0.5,-0.5,-0.2,-0.1, 0.0,-0.4,-0.1, 0.7,-0.4,-0.3,-0.3, 0.0, 0.3,-0.2,-0.4,-0.4,-0.2,-0.6}, //PHE
        {-0.1,-0.3,-0.3,-0.3,-0.4,-0.2,-0.2,-0.3,-0.3,-0.4,-0.4,-0.2,-0.3,-0.4, 0.8,-0.2,-0.2,-0.5,-0.4,-0.3,-0.3,-0.2,-0.2,-0.6}, //PRO
        { 0.1,-0.1, 0.0,-0.1,-0.2,-0.1,-0.1,-0.1,-0.2,-0.3,-0.3,-0.1,-0.2,-0.3,-0.2, 0.5, 0.1,-0.4,-0.3,-0.2, 0.0,-0.1,-0.1,-0.6}, //SER
        { 0.0,-0.2, 0.0,-0.2,-0.2,-0.1,-0.1,-0.3,-0.2,-0.1,-0.2,-0.1,-0.1,-0.3,-0.2, 0.1, 0.6,-0.4,-0.2,-0.1,-0.1,-0.1,-0.1,-0.6}, //THR
        {-0.4,-0.4,-0.5,-0.6,-0.4,-0.3,-0.5,-0.4,-0.3,-0.4,-0.3,-0.5,-0.2, 0.0,-0.5,-0.4,-0.4, 1.1, 0.2,-0.3,-0.6,-0.4,-0.3,-0.6}, //TRP
        {-0.3,-0.3,-0.3,-0.4,-0.4,-0.3,-0.4,-0.5, 0.1,-0.2,-0.2,-0.3,-0.2, 0.3,-0.4,-0.3,-0.2, 0.2, 0.8,-0.3,-0.4,-0.3,-0.2,-0.6}, //TYR
        {-0.1,-0.3,-0.4,-0.5,-0.2,-0.3,-0.3,-0.5,-0.4, 0.3, 0.0,-0.3, 0.0,-0.2,-0.3,-0.2,-0.1,-0.3,-0.3, 0.5,-0.4,-0.3,-0.2,-0.6}, //VAL
        {-0.2,-0.2, 0.4, 0.4,-0.4,-0.1, 0.0,-0.2,-0.1,-0.5,-0.5,-0.1,-0.4,-0.4,-0.3, 0.0,-0.1,-0.6,-0.4,-0.4, 0.4, 0.0,-0.2,-0.6}, //ASX
        {-0.1, 0.0,-0.1, 0.0,-0.5, 0.4, 0.4,-0.3, 0.0,-0.4,-0.4, 0.1,-0.2,-0.4,-0.2,-0.1,-0.1,-0.4,-0.3,-0.3, 0.0, 0.4,-0.1,-0.6}, //GLX
        {-0.1,-0.2,-0.2,-0.2,-0.3,-0.1,-0.2,-0.2,-0.2,-0.2,-0.2,-0.1,-0.1,-0.2,-0.2,-0.1,-0.1,-0.3,-0.2,-0.2,-0.2,-0.1,-0.2,-0.6}, //XXX
        {-0.6,-0.6,-0.6,-0.6,-0.6,-0.6,-0.6,-0.6,-0.6,-0.6,-0.6,-0.6,-0.6,-0.6,-0.6,-0.6,-0.6,-0.6,-0.6,-0.6,-0.6,-0.6,-0.6, 0.1}  //GAP
      },
      //BLOSUM_80
      {
//        ALA  ARG  ASN  ASP  CYS  GLN  GLU  GLY  HIS  ILE  LEU  LYS  MET  PHE  PRO  SER  THR  TRP  TYR  VAL  ASX  GLX  XXX  GAP
        { 0.7,-0.3,-0.3,-0.3,-0.1,-0.2,-0.2, 0.0,-0.3,-0.3,-0.3,-0.1,-0.2,-0.4,-0.1, 0.2, 0.0,-0.5,-0.4,-0.1,-0.3,-0.2,-0.1,-0.8}, //ALA
        {-0.3, 0.9,-0.1,-0.3,-0.6, 0.1,-0.1,-0.4, 0.0,-0.5,-0.4, 0.3,-0.3,-0.5,-0.3,-0.2,-0.2,-0.5,-0.4,-0.4,-0.2, 0.0,-0.2,-0.8}, //ARG
        {-0.3,-0.1, 0.9, 0.2,-0.5, 0.0,-0.1,-0.1, 0.1,-0.6,-0.6, 0.0,-0.4,-0.6,-0.4, 0.1, 0.0,-0.7,-0.4,-0.5, 0.5,-0.1,-0.2,-0.8}, //ASN
        {-0.3,-0.3, 0.2, 1.0,-0.7,-0.1, 0.2,-0.3,-0.2,-0.7,-0.7,-0.2,-0.6,-0.6,-0.3,-0.1,-0.2,-0.8,-0.6,-0.6, 0.6, 0.1,-0.3,-0.8}, //ASP
        {-0.1,-0.6,-0.5,-0.7, 1.3,-0.5,-0.7,-0.6,-0.7,-0.2,-0.3,-0.6,-0.3,-0.4,-0.6,-0.2,-0.2,-0.5,-0.5,-0.2,-0.6,-0.7,-0.4,-0.8}, //CYS
        {-0.2, 0.1, 0.0,-0.1,-0.5, 0.9, 0.3,-0.4, 0.1,-0.5,-0.4, 0.2,-0.1,-0.5,-0.3,-0.1,-0.1,-0.4,-0.3,-0.4,-0.1, 0.5,-0.2,-0.8}, //GLN
        {-0.2,-0.1,-0.1, 0.2,-0.7, 0.3, 0.8,-0.4, 0.0,-0.6,-0.6, 0.1,-0.4,-0.6,-0.2,-0.1,-0.2,-0.6,-0.5,-0.4, 0.1, 0.6,-0.2,-0.8}, //GLU
        { 0.0,-0.4,-0.1,-0.3,-0.6,-0.4,-0.4, 0.9,-0.4,-0.7,-0.7,-0.3,-0.5,-0.6,-0.5,-0.1,-0.3,-0.6,-0.6,-0.6,-0.2,-0.4,-0.3,-0.8}, //GLY
        {-0.3, 0.0, 0.1,-0.2,-0.7, 0.1, 0.0,-0.4, 1.2,-0.6,-0.5,-0.1,-0.4,-0.2,-0.4,-0.2,-0.3,-0.4, 0.3,-0.5,-0.1, 0.0,-0.2,-0.8}, //HIS
        {-0.3,-0.5,-0.6,-0.7,-0.2,-0.5,-0.6,-0.7,-0.6, 0.7, 0.2,-0.5, 0.2,-0.1,-0.5,-0.4,-0.2,-0.5,-0.3, 0.4,-0.6,-0.6,-0.2,-0.8}, //ILE
        {-0.3,-0.4,-0.6,-0.7,-0.3,-0.4,-0.6,-0.7,-0.5, 0.2, 0.6,-0.4, 0.3, 0.0,-0.5,-0.4,-0.3,-0.4,-0.2, 0.1,-0.7,-0.5,-0.2,-0.8}, //LEU
        {-0.1, 0.3, 0.0,-0.2,-0.6, 0.2, 0.1,-0.3,-0.1,-0.5,-0.4, 0.8,-0.3,-0.5,-0.2,-0.1,-0.1,-0.6,-0.4,-0.4,-0.1, 0.1,-0.2,-0.8}, //LYS
        {-0.2,-0.3,-0.4,-0.6,-0.3,-0.1,-0.4,-0.5,-0.4, 0.2, 0.3,-0.3, 0.9, 0.0,-0.4,-0.3,-0.1,-0.3,-0.3, 0.1,-0.5,-0.3,-0.2,-0.8}, //MET
        {-0.4,-0.5,-0.6,-0.6,-0.4,-0.5,-0.6,-0.6,-0.2,-0.1, 0.0,-0.5, 0.0, 1.0,-0.6,-0.4,-0.4, 0.0, 0.4,-0.2,-0.6,-0.6,-0.3,-0.8}, //PHE
        {-0.1,-0.3,-0.4,-0.3,-0.6,-0.3,-0.2,-0.5,-0.4,-0.5,-0.5,-0.2,-0.4,-0.6, 1.2,-0.2,-0.3,-0.7,-0.6,-0.4,-0.4,-0.2,-0.3,-0.8}, //PRO
        { 0.2,-0.2, 0.1,-0.1,-0.2,-0.1,-0.1,-0.1,-0.2,-0.4,-0.4,-0.1,-0.3,-0.4,-0.2, 0.7, 0.2,-0.6,-0.3,-0.3, 0.0,-0.1,-0.1,-0.8}, //SER
        { 0.0,-0.2, 0.0,-0.2,-0.2,-0.1,-0.2,-0.3,-0.3,-0.2,-0.3,-0.1,-0.1,-0.4,-0.3, 0.2, 0.8,-0.5,-0.3, 0.0,-0.1,-0.2,-0.1,-0.8}, //THR
        {-0.5,-0.5,-0.7,-0.8,-0.5,-0.4,-0.6,-0.6,-0.4,-0.5,-0.4,-0.6,-0.3, 0.0,-0.7,-0.6,-0.5, 1.6, 0.3,-0.5,-0.8,-0.5,-0.5,-0.8}, //TRP
        {-0.4,-0.4,-0.4,-0.6,-0.5,-0.3,-0.5,-0.6, 0.3,-0.3,-0.2,-0.4,-0.3, 0.4,-0.6,-0.3,-0.3, 0.3, 1.1,-0.3,-0.5,-0.4,-0.3,-0.8}, //TYR
        {-0.1,-0.4,-0.5,-0.6,-0.2,-0.4,-0.4,-0.6,-0.5, 0.4, 0.1,-0.4, 0.1,-0.2,-0.4,-0.3, 0.0,-0.5,-0.3, 0.7,-0.6,-0.4,-0.2,-0.8}, //VAL
        {-0.3,-0.2, 0.5, 0.6,-0.6,-0.1, 0.1,-0.2,-0.1,-0.6,-0.7,-0.1,-0.5,-0.6,-0.4, 0.0,-0.1,-0.8,-0.5,-0.6, 0.6, 0.0,-0.3,-0.8}, //ASX
        {-0.2, 0.0,-0.1, 0.1,-0.7, 0.5, 0.6,-0.4, 0.0,-0.6,-0.5, 0.1,-0.3,-0.6,-0.2,-0.1,-0.2,-0.5,-0.4,-0.4, 0.0, 0.6,-0.1,-0.8}, //GLX
        {-0.1,-0.2,-0.2,-0.3,-0.4,-0.2,-0.2,-0.3,-0.2,-0.2,-0.2,-0.2,-0.2,-0.3,-0.3,-0.1,-0.1,-0.5,-0.3,-0.2,-0.3,-0.1,-0.2,-0.8}, //XXX
        {-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8, 0.1}  //GAP
      },
      //BLOSUM_62
      {
//        ALA  ARG  ASN  ASP  CYS  GLN  GLU  GLY  HIS  ILE  LEU  LYS  MET  PHE  PRO  SER  THR  TRP  TYR  VAL  ASX  GLX  XXX  GAP
        { 0.4,-0.1,-0.2,-0.2, 0.0,-0.1,-0.1, 0.0,-0.2,-0.1,-0.1,-0.1,-0.1,-0.2,-0.1, 0.1, 0.0,-0.3,-0.2, 0.0,-0.2,-0.1, 0.0,-0.4}, //ALA
        {-0.1, 0.5, 0.0,-0.2,-0.3, 0.1, 0.0,-0.2, 0.0,-0.3,-0.2, 0.2,-0.1,-0.3,-0.2,-0.1,-0.1,-0.3,-0.2,-0.3,-0.1, 0.0,-0.1,-0.4}, //ARG
        {-0.2, 0.0, 0.6, 0.1,-0.3, 0.0, 0.0, 0.0, 0.1,-0.3,-0.3, 0.0,-0.2,-0.3,-0.2, 0.1, 0.0,-0.4,-0.2,-0.3, 0.3, 0.0,-0.1,-0.4}, //ASN
        {-0.2,-0.2, 0.1, 0.6,-0.3, 0.0, 0.2,-0.1,-0.1,-0.3,-0.4,-0.1,-0.3,-0.3,-0.1, 0.0,-0.1,-0.4,-0.3,-0.3, 0.4, 0.1,-0.1,-0.4}, //ASP
        { 0.0,-0.3,-0.3,-0.3, 0.9,-0.3,-0.4,-0.3,-0.3,-0.1,-0.1,-0.3,-0.1,-0.2,-0.3,-0.1,-0.1,-0.2,-0.2,-0.1,-0.3,-0.3,-0.2,-0.4}, //CYS
        {-0.1, 0.1, 0.0, 0.0,-0.3, 0.5, 0.2,-0.2, 0.0,-0.3,-0.2, 0.1, 0.0,-0.3,-0.1, 0.0,-0.1,-0.2,-0.1,-0.2, 0.0, 0.3,-0.1,-0.4}, //GLN
        {-0.1, 0.0, 0.0, 0.2,-0.4, 0.2, 0.5,-0.2, 0.0,-0.3,-0.3, 0.1,-0.2,-0.3,-0.1, 0.0,-0.1,-0.3,-0.2,-0.2, 0.1, 0.4,-0.1,-0.4}, //GLU
        { 0.0,-0.2, 0.0,-0.1,-0.3,-0.2,-0.2, 0.6,-0.2,-0.4,-0.4,-0.2,-0.3,-0.3,-0.2, 0.0,-0.2,-0.2,-0.3,-0.3,-0.1,-0.2,-0.1,-0.4}, //GLY
        {-0.2, 0.0, 0.1,-0.1,-0.3, 0.0, 0.0,-0.2, 0.8,-0.3,-0.3,-0.1,-0.2,-0.1,-0.2,-0.1,-0.2,-0.2, 0.2,-0.3, 0.0, 0.0,-0.1,-0.4}, //HIS
        {-0.1,-0.3,-0.3,-0.3,-0.1,-0.3,-0.3,-0.4,-0.3, 0.4, 0.2,-0.3, 0.1, 0.0,-0.3,-0.2,-0.1,-0.3,-0.1, 0.3,-0.3,-0.3,-0.1,-0.4}, //ILE
        {-0.1,-0.2,-0.3,-0.4,-0.1,-0.2,-0.3,-0.4,-0.3, 0.2, 0.4,-0.2, 0.2, 0.0,-0.3,-0.2,-0.1,-0.2,-0.1, 0.1,-0.4,-0.3,-0.1,-0.4}, //LEU
        {-0.1, 0.2, 0.0,-0.1,-0.3, 0.1, 0.1,-0.2,-0.1,-0.3,-0.2, 0.5,-0.1,-0.3,-0.1, 0.0,-0.1,-0.3,-0.2,-0.2, 0.0, 0.1,-0.1,-0.4}, //LYS
        {-0.1,-0.1,-0.2,-0.3,-0.1, 0.0,-0.2,-0.3,-0.2, 0.1, 0.2,-0.1, 0.5, 0.0,-0.2,-0.1,-0.1,-0.1,-0.1, 0.1,-0.3,-0.1,-0.1,-0.4}, //MET
        {-0.2,-0.3,-0.3,-0.3,-0.2,-0.3,-0.3,-0.3,-0.1, 0.0, 0.0,-0.3, 0.0, 0.6,-0.4,-0.2,-0.2, 0.1, 0.3,-0.1,-0.3,-0.3,-0.1,-0.4}, //PHE
        {-0.1,-0.2,-0.2,-0.1,-0.3,-0.1,-0.1,-0.2,-0.2,-0.3,-0.3,-0.1,-0.2,-0.4, 0.7,-0.1,-0.1,-0.4,-0.3,-0.2,-0.2,-0.1,-0.2,-0.4}, //PRO
        { 0.1,-0.1, 0.1, 0.0,-0.1, 0.0, 0.0, 0.0,-0.1,-0.2,-0.2, 0.0,-0.1,-0.2,-0.1, 0.4, 0.1,-0.3,-0.2,-0.2, 0.0, 0.0, 0.0,-0.4}, //SER
        { 0.0,-0.1, 0.0,-0.1,-0.1,-0.1,-0.1,-0.2,-0.2,-0.1,-0.1,-0.1,-0.1,-0.2,-0.1, 0.1, 0.5,-0.2,-0.2, 0.0,-0.1,-0.1, 0.0,-0.4}, //THR
        {-0.3,-0.3,-0.4,-0.4,-0.2,-0.2,-0.3,-0.2,-0.2,-0.3,-0.2,-0.3,-0.1, 0.1,-0.4,-0.3,-0.2, 1.1, 0.2,-0.3,-0.4,-0.3,-0.2,-0.4}, //TRP
        {-0.2,-0.2,-0.2,-0.3,-0.2,-0.1,-0.2,-0.3, 0.2,-0.1,-0.1,-0.2,-0.1, 0.3,-0.3,-0.2,-0.2, 0.2, 0.7,-0.1,-0.3,-0.2,-0.1,-0.4}, //TYR
        { 0.0,-0.3,-0.3,-0.3,-0.1,-0.2,-0.2,-0.3,-0.3, 0.3, 0.1,-0.2, 0.1,-0.1,-0.2,-0.2, 0.0,-0.3,-0.1, 0.4,-0.3,-0.2,-0.1,-0.4}, //VAL
        {-0.2,-0.1, 0.3, 0.4,-0.3, 0.0, 0.1,-0.1, 0.0,-0.3,-0.4, 0.0,-0.3,-0.3,-0.2, 0.0,-0.1,-0.4,-0.3,-0.3, 0.4, 0.1,-0.1,-0.4}, //ASX
        {-0.1, 0.0, 0.0, 0.1,-0.3, 0.3, 0.4,-0.2, 0.0,-0.3,-0.3, 0.1,-0.1,-0.3,-0.1, 0.0,-0.1,-0.3,-0.2,-0.2, 0.1, 0.4,-0.1,-0.4}, //GLX
        { 0.0,-0.1,-0.1,-0.1,-0.2,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.2, 0.0, 0.0,-0.2,-0.1,-0.1,-0.1,-0.1,-0.1,-0.4}, //XXX
        {-0.4,-0.4,-0.4,-0.4,-0.4,-0.4,-0.4,-0.4,-0.4,-0.4,-0.4,-0.4,-0.4,-0.4,-0.4,-0.4,-0.4,-0.4,-0.4,-0.4,-0.4,-0.4,-0.4, 0.1}  //GAP
      },
      //BLOSUM_45
      {
//        ALA  ARG  ASN  ASP  CYS  GLN  GLU  GLY  HIS  ILE  LEU  LYS  MET  PHE  PRO  SER  THR  TRP  TYR  VAL  ASX  GLX  XXX  GAP
        { 0.5,-0.2,-0.1,-0.2,-0.1,-0.1,-0.1, 0.0,-0.2,-0.1,-0.1,-0.1,-0.1,-0.2,-0.1, 0.1, 0.0,-0.2,-0.2, 0.0,-0.1,-0.1, 0.0,-0.5}, //ALA
        {-0.2, 0.7, 0.0,-0.1,-0.3, 0.1, 0.0,-0.2, 0.0,-0.3,-0.2, 0.3,-0.1,-0.2,-0.2,-0.1,-0.1,-0.2,-0.1,-0.2,-0.1, 0.0,-0.1,-0.5}, //ARG
        {-0.1, 0.0, 0.6, 0.2,-0.2, 0.0, 0.0, 0.0, 0.1,-0.2,-0.3, 0.0,-0.2,-0.2,-0.2, 0.1, 0.0,-0.4,-0.2,-0.3, 0.4, 0.0,-0.1,-0.5}, //ASN
        {-0.2,-0.1, 0.2, 0.7,-0.3, 0.0, 0.2,-0.1, 0.0,-0.4,-0.3, 0.0,-0.3,-0.4,-0.1, 0.0,-0.1,-0.4,-0.2,-0.3, 0.5, 0.1,-0.1,-0.5}, //ASP
        {-0.1,-0.3,-0.2,-0.3, 1.2,-0.3,-0.3,-0.3,-0.3,-0.3,-0.2,-0.3,-0.2,-0.2,-0.4,-0.1,-0.1,-0.5,-0.3,-0.1,-0.2,-0.3,-0.2,-0.5}, //CYS
        {-0.1, 0.1, 0.0, 0.0,-0.3, 0.6, 0.2,-0.2, 0.1,-0.2,-0.2, 0.1, 0.0,-0.4,-0.1, 0.0,-0.1,-0.2,-0.1,-0.3, 0.0, 0.4,-0.1,-0.5}, //GLN
        {-0.1, 0.0, 0.0, 0.2,-0.3, 0.2, 0.6,-0.2, 0.0,-0.3,-0.2, 0.1,-0.2,-0.3, 0.0, 0.0,-0.1,-0.3,-0.2,-0.3, 0.1, 0.4,-0.1,-0.5}, //GLU
        { 0.0,-0.2, 0.0,-0.1,-0.3,-0.2,-0.2, 0.7,-0.2,-0.4,-0.3,-0.2,-0.2,-0.3,-0.2, 0.0,-0.2,-0.2,-0.3,-0.3,-0.1,-0.2,-0.1,-0.5}, //GLY
        {-0.2, 0.0, 0.1, 0.0,-0.3, 0.1, 0.0,-0.2, 1.0,-0.3,-0.2,-0.1, 0.0,-0.2,-0.2,-0.1,-0.2,-0.3, 0.2,-0.3, 0.0, 0.0,-0.1,-0.5}, //HIS
        {-0.1,-0.3,-0.2,-0.4,-0.3,-0.2,-0.3,-0.4,-0.3, 0.5, 0.2,-0.3, 0.2, 0.0,-0.2,-0.2,-0.1,-0.2, 0.0, 0.3,-0.3,-0.3,-0.1,-0.5}, //ILE
        {-0.1,-0.2,-0.3,-0.3,-0.2,-0.2,-0.2,-0.3,-0.2, 0.2, 0.5,-0.3, 0.2, 0.1,-0.3,-0.3,-0.1,-0.2, 0.0, 0.1,-0.3,-0.2,-0.1,-0.5}, //LEU
        {-0.1, 0.3, 0.0, 0.0,-0.3, 0.1, 0.1,-0.2,-0.1,-0.3,-0.3, 0.5,-0.1,-0.3,-0.1,-0.1,-0.1,-0.2,-0.1,-0.2, 0.0, 0.1,-0.1,-0.5}, //LYS
        {-0.1,-0.1,-0.2,-0.3,-0.2, 0.0,-0.2,-0.2, 0.0, 0.2, 0.2,-0.1, 0.6, 0.0,-0.2,-0.2,-0.1,-0.2, 0.0, 0.1,-0.2,-0.1,-0.1,-0.5}, //MET
        {-0.2,-0.2,-0.2,-0.4,-0.2,-0.4,-0.3,-0.3,-0.2, 0.0, 0.1,-0.3, 0.0, 0.8,-0.3,-0.2,-0.1, 0.1, 0.3, 0.0,-0.3,-0.3,-0.1,-0.5}, //PHE
        {-0.1,-0.2,-0.2,-0.1,-0.4,-0.1, 0.0,-0.2,-0.2,-0.2,-0.3,-0.1,-0.2,-0.3, 0.9,-0.1,-0.1,-0.3,-0.3,-0.3,-0.2,-0.1,-0.1,-0.5}, //PRO
        { 0.1,-0.1, 0.1, 0.0,-0.1, 0.0, 0.0, 0.0,-0.1,-0.2,-0.3,-0.1,-0.2,-0.2,-0.1, 0.4, 0.2,-0.4,-0.2,-0.1, 0.0, 0.0, 0.0,-0.5}, //SER
        { 0.0,-0.1, 0.0,-0.1,-0.1,-0.1,-0.1,-0.2,-0.2,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1, 0.2, 0.5,-0.3,-0.1, 0.0, 0.0,-0.1, 0.0,-0.5}, //THR
        {-0.2,-0.2,-0.4,-0.4,-0.5,-0.2,-0.3,-0.2,-0.3,-0.2,-0.2,-0.2,-0.2, 0.1,-0.3,-0.4,-0.3, 1.5, 0.3,-0.3,-0.4,-0.2,-0.2,-0.5}, //TRP
        {-0.2,-0.1,-0.2,-0.2,-0.3,-0.1,-0.2,-0.3, 0.2, 0.0, 0.0,-0.1, 0.0, 0.3,-0.3,-0.2,-0.1, 0.3, 0.8,-0.1,-0.2,-0.2,-0.1,-0.5}, //TYR
        { 0.0,-0.2,-0.3,-0.3,-0.1,-0.3,-0.3,-0.3,-0.3, 0.3, 0.1,-0.2, 0.1, 0.0,-0.3,-0.1, 0.0,-0.3,-0.1, 0.5,-0.3,-0.3,-0.1,-0.5}, //VAL
        {-0.1,-0.1, 0.4, 0.5,-0.2, 0.0, 0.1,-0.1, 0.0,-0.3,-0.3, 0.0,-0.2,-0.3,-0.2, 0.0, 0.0,-0.4,-0.2,-0.3, 0.4, 0.2,-0.1,-0.5}, //ASX
        {-0.1, 0.0, 0.0, 0.1,-0.3, 0.4, 0.4,-0.2, 0.0,-0.3,-0.2, 0.1,-0.1,-0.3,-0.1, 0.0,-0.1,-0.2,-0.2,-0.3, 0.2, 0.4,-0.1,-0.5}, //GLX
        { 0.0,-0.1,-0.1,-0.1,-0.2,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1, 0.0, 0.0,-0.2,-0.1,-0.1,-0.1,-0.1,-0.1,-0.5}, //XXX
        {-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5, 0.1}  //GAP
      }
    };

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a BLOSUM table
    //! @param BLOSUM_TABLE TableType to be used
    AAAssignmentBLOSUM::AAAssignmentBLOSUM( const TableType &BLOSUM_TABLE) :
      m_BLOSUMTable( BLOSUM_TABLE)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to a new AAAssignmentBLOSUM copied from this one
    AAAssignmentBLOSUM *AAAssignmentBLOSUM::Clone() const
    {
      return new AAAssignmentBLOSUM( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AAAssignmentBLOSUM::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief return the mutation probability for two AAs given a table
    //! @param BLOSUM_TABLE TableType to be used
    //! @param AA_TYPE_A first AAType of interest
    //! @param AA_TYPE_B second AAType of interest
    //! @return the mutation probability for two AAs given a table
    double AAAssignmentBLOSUM::Probability
    (
      const TableType BLOSUM_TABLE,
      const biol::AAType &AA_TYPE_A,
      const biol::AAType &AA_TYPE_B
    )
    {
      return pow( 10, s_BLOSUMTable[ BLOSUM_TABLE][ AA_TYPE_A][ AA_TYPE_B]) / double( 20);
    }

    //! @brief returns the requested BLOSUM table as a matrix
    //! @param BLOSUM_TABLE requested BLOSUM table name
    //! @return the requested BLOSUM table as a matrix
    linal::Matrix< double> AAAssignmentBLOSUM::GetBLOSUMMatrix( const TableType &BLOSUM_TABLE)
    {
      // construct the matrix and return it
      return linal::Matrix< double>
      (
        biol::AATypes::s_NumberStandardAATypes + 4,
        biol::AATypes::s_NumberStandardAATypes + 4,
        &s_BLOSUMTable[ BLOSUM_TABLE][ 0][ 0]
      );
    }

    //! @brief returns the requested BLOSUM table as a matrix
    //! @param BLOSUM_TABLE requested BLOSUM table name
    //! @return the requested BLOSUM table as a matrix
    linal::VectorConstReference< double> AAAssignmentBLOSUM::GetBLOSUMRow
    (
      const TableType &BLOSUM_TABLE,
      const biol::AAType &AATYPE
    )
    {
      return
        linal::VectorConstReference< double>
        (
          biol::AATypes::s_NumberStandardAATypes + 4,
          s_BLOSUMTable[ BLOSUM_TABLE][ AATYPE]
        );
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that calculates the score between two assigned members
    //! @param MEMBER_A amino acid A that is compared
    //! @param MEMBER_B amino acid A that is compared
    //! @return value from the m_BLOSUMTable for this combination of amino acids
    double AAAssignmentBLOSUM::operator()( const biol::AABase &MEMBER_A, const biol::AABase &MEMBER_B) const
    {
      static const size_t matrix_dim( biol::AATypes::s_NumberStandardAATypes + 4);

      // check that the types are valid, if not, return undefined
      if( MEMBER_A.GetType() >= matrix_dim || MEMBER_B.GetType() >= matrix_dim)
      {
        return util::GetUndefined< double>();
      }

      return s_BLOSUMTable[ m_BLOSUMTable][ MEMBER_A.GetType()][ MEMBER_B.GetType()];
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write to ostream
    //! @param OSTREAM is the output stream
    //! @param INDENT indentation
    //! @return returns the output stream
    std::ostream &AAAssignmentBLOSUM::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_BLOSUMTable, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief read from istream
    //! @param ISTREAM is the input stream
    //! @return returns the input stream
    std::istream &AAAssignmentBLOSUM::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_BLOSUMTable, ISTREAM);

      // end
      return ISTREAM;
    }

  } // namespace score
} // namespace bcl

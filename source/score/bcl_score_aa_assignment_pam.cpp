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
#include "score/bcl_score_aa_assignment_pam.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  ///////////
  // enums //
  ///////////

    //! @brief TableType as string
    //! @param TABLE_TYPE the TableType
    //! @return the string for the TableType
    const std::string &AAAssignmentPAM::GetTableDescriptor( const TableType &TABLE_TYPE)
    {
      static const std::string s_descriptors[] =
      {
        "PAM_100",
        "PAM_120",
        "PAM_160",
        "PAM_250",
        GetStaticClassName< TableType>()
      };

      return s_descriptors[ TABLE_TYPE];
    }

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AAAssignmentPAM::s_Instance
    (
      GetObjectInstances().AddInstance( new AAAssignmentPAM())
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
    //! PAM (Point Accepted Mutation) matrix:
    //! Amino acid scoring matrices are traditionally PAM (Point Accepted Mutation)
    //! matrices which refer to various degrees of sensitivity depending on the evolutionary distance
    //! between sequence pairs. In this manner PAM40 is most sensitive for sequences 40 PAMs
    //! apart. PAM250 is for more distantly related sequences and is considered a good general matrix
    //! for protein database searching. For nucleotide sequence searching a simpler approach is used which either
    //! convert a PAM40 matrix into match/mismatch values which takes into consideration that a purine may be
    //! replaced by a purine and a pyrimidine by a pyrimidine.
    const double AAAssignmentPAM::s_PAMTable[][ biol::AATypes::s_NumberStandardAATypes + 4][ biol::AATypes::s_NumberStandardAATypes + 4] =
    {
      //PAM_100
      {
//        ALA  ARG  ASN  ASP  CYS  GLN  GLU  GLY  HIS  ILE  LEU  LYS  MET  PHE  PRO  SER  THR  TRP  TYR  VAL  ASX  GLX  XXX  GAP
        { 0.4,-0.3,-0.1,-0.1,-0.3,-0.2, 0.0, 0.1,-0.3,-0.2,-0.3,-0.3,-0.2,-0.5, 0.1, 0.1, 0.1,-0.7,-0.4, 0.0,-0.1,-0.1,-0.1,-0.9},//ALA
        {-0.3, 0.7,-0.2,-0.4,-0.5, 0.1,-0.3,-0.5, 0.1,-0.3,-0.5, 0.2,-0.1,-0.6,-0.1,-0.1,-0.3, 0.1,-0.6,-0.4,-0.3,-0.1,-0.2,-0.9},//ARG
        {-0.1,-0.2, 0.5, 0.3,-0.5,-0.1, 0.1,-0.1, 0.2,-0.3,-0.4, 0.1,-0.4,-0.5,-0.2, 0.1, 0.0,-0.5,-0.2,-0.3, 0.4, 0.0,-0.1,-0.9},//ASN
        {-0.1,-0.4, 0.3, 0.5,-0.7, 0.0, 0.4,-0.1,-0.1,-0.4,-0.6,-0.1,-0.5,-0.8,-0.3,-0.1,-0.2,-0.9,-0.6,-0.4, 0.4, 0.3,-0.2,-0.9},//ASP
        {-0.3,-0.5,-0.5,-0.7, 0.9,-0.8,-0.8,-0.5,-0.4,-0.3,-0.8,-0.8,-0.7,-0.7,-0.4,-0.1,-0.4,-0.9,-0.1,-0.3,-0.6,-0.8,-0.5,-0.9},//CYS
        {-0.2, 0.1,-0.1, 0.0,-0.8, 0.6, 0.2,-0.3, 0.3,-0.4,-0.2, 0.0,-0.2,-0.7,-0.1,-0.2,-0.2,-0.7,-0.6,-0.3, 0.0, 0.5,-0.2,-0.9},//GLN
        { 0.0,-0.3, 0.1, 0.4,-0.8, 0.2, 0.5,-0.1,-0.1,-0.3,-0.5,-0.1,-0.4,-0.8,-0.2,-0.1,-0.2,-0.9,-0.5,-0.3, 0.3, 0.4,-0.2,-0.9},//GLU
        { 0.1,-0.5,-0.1,-0.1,-0.5,-0.3,-0.1, 0.5,-0.4,-0.5,-0.6,-0.3,-0.4,-0.6,-0.2, 0.0,-0.2,-0.9,-0.7,-0.3,-0.1,-0.2,-0.2,-0.9},//GLY
        {-0.3, 0.1, 0.2,-0.1,-0.4, 0.3,-0.1,-0.4, 0.7,-0.4,-0.3,-0.2,-0.4,-0.3,-0.1,-0.2,-0.3,-0.4,-0.1,-0.3, 0.1, 0.1,-0.2,-0.9},//HIS
        {-0.2,-0.3,-0.3,-0.4,-0.3,-0.4,-0.3,-0.5,-0.4, 0.6, 0.1,-0.3, 0.1, 0.0,-0.4,-0.3, 0.0,-0.7,-0.3, 0.3,-0.3,-0.3,-0.2,-0.9},//ILE
        {-0.3,-0.5,-0.4,-0.6,-0.8,-0.2,-0.5,-0.6,-0.3, 0.1, 0.6,-0.4, 0.3, 0.0,-0.4,-0.4,-0.3,-0.3,-0.3, 0.0,-0.5,-0.4,-0.3,-0.9},//LEU
        {-0.3, 0.2, 0.1,-0.1,-0.8, 0.0,-0.1,-0.3,-0.2,-0.3,-0.4, 0.5, 0.0,-0.7,-0.3,-0.1,-0.1,-0.6,-0.6,-0.4, 0.0,-0.1,-0.2,-0.9},//LYS
        {-0.2,-0.1,-0.4,-0.5,-0.7,-0.2,-0.4,-0.4,-0.4, 0.1, 0.3, 0.0, 0.9,-0.1,-0.4,-0.3,-0.1,-0.6,-0.5, 0.1,-0.4,-0.2,-0.2,-0.9},//MET
        {-0.5,-0.6,-0.5,-0.8,-0.7,-0.7,-0.8,-0.6,-0.3, 0.0, 0.0,-0.7,-0.1, 0.8,-0.6,-0.4,-0.5,-0.1, 0.4,-0.3,-0.6,-0.7,-0.4,-0.9},//PHE
        { 0.1,-0.1,-0.2,-0.3,-0.4,-0.1,-0.2,-0.2,-0.1,-0.4,-0.4,-0.3,-0.4,-0.6, 0.7, 0.0,-0.1,-0.7,-0.7,-0.3,-0.3,-0.1,-0.2,-0.9},//PRO
        { 0.1,-0.1, 0.1,-0.1,-0.1,-0.2,-0.1, 0.0,-0.2,-0.3,-0.4,-0.1,-0.3,-0.4, 0.0, 0.4, 0.2,-0.3,-0.4,-0.2, 0.0,-0.2,-0.1,-0.9},//SER
        { 0.1,-0.3, 0.0,-0.2,-0.4,-0.2,-0.2,-0.2,-0.3, 0.0,-0.3,-0.1,-0.1,-0.5,-0.1, 0.2, 0.5,-0.7,-0.4, 0.0,-0.1,-0.2,-0.1,-0.9},//THR
        {-0.7, 0.1,-0.5,-0.9,-0.9,-0.7,-0.9,-0.9,-0.4,-0.7,-0.3,-0.6,-0.6,-0.1,-0.7,-0.3,-0.7, 1.2,-0.2,-0.9,-0.6,-0.8,-0.6,-0.9},//TRP
        {-0.4,-0.6,-0.2,-0.6,-0.1,-0.6,-0.5,-0.7,-0.1,-0.3,-0.3,-0.6,-0.5, 0.4,-0.7,-0.4,-0.4,-0.2, 0.9,-0.4,-0.4,-0.6,-0.4,-0.9},//TYR
        { 0.0,-0.4,-0.3,-0.4,-0.3,-0.3,-0.3,-0.3,-0.3, 0.3, 0.0,-0.4, 0.1,-0.3,-0.3,-0.2, 0.0,-0.9,-0.4, 0.5,-0.4,-0.3,-0.2,-0.9},//VAL
        {-0.1,-0.3, 0.4, 0.4,-0.6, 0.0, 0.3,-0.1, 0.1,-0.3,-0.5, 0.0,-0.4,-0.6,-0.3, 0.0,-0.1,-0.6,-0.4,-0.4, 0.4, 0.2,-0.2,-0.9},//ASX
        {-0.1,-0.1, 0.0, 0.3,-0.8, 0.5, 0.4,-0.2, 0.1,-0.3,-0.4,-0.1,-0.2,-0.7,-0.1,-0.2,-0.2,-0.8,-0.6,-0.3, 0.2, 0.5,-0.2,-0.9},//GLX
        {-0.1,-0.2,-0.1,-0.2,-0.5,-0.2,-0.2,-0.2,-0.2,-0.2,-0.3,-0.2,-0.2,-0.4,-0.2,-0.1,-0.1,-0.6,-0.4,-0.2,-0.2,-0.2,-0.2,-0.9},//XXX
        {-0.9,-0.9,-0.9,-0.9,-0.9,-0.9,-0.9,-0.9,-0.9,-0.9,-0.9,-0.9,-0.9,-0.9,-0.9,-0.9,-0.9,-0.9,-0.9,-0.9,-0.9,-0.9,-0.9, 0.1} //GAP
      },
      //PAM_120
      {
//        ALA  ARG  ASN  ASP  CYS  GLN  GLU  GLY  HIS  ILE  LEU  LYS  MET  PHE  PRO  SER  THR  TRP  TYR  VAL  ASX  GLX  XXX  GAP
        { 0.3,-0.3,-0.1, 0.0,-0.3,-0.1, 0.0, 0.1,-0.3,-0.1,-0.3,-0.2,-0.2,-0.4, 0.1, 0.1, 0.1,-0.7,-0.4, 0.0, 0.0,-0.1,-0.1,-0.8},//ALA
        {-0.3, 0.6,-0.1,-0.3,-0.4, 0.1,-0.3,-0.4, 0.1,-0.2,-0.4, 0.2,-0.1,-0.5,-0.1,-0.1,-0.2, 0.1,-0.5,-0.3,-0.2,-0.1,-0.2,-0.8},//ARG
        {-0.1,-0.1, 0.4, 0.2,-0.5, 0.0, 0.1, 0.0, 0.2,-0.2,-0.4, 0.1,-0.3,-0.4,-0.2, 0.1, 0.0,-0.4,-0.2,-0.3, 0.3, 0.0,-0.1,-0.8},//ASN
        { 0.0,-0.3, 0.2, 0.5,-0.7, 0.1, 0.3, 0.0, 0.0,-0.3,-0.5,-0.1,-0.4,-0.7,-0.3, 0.0,-0.1,-0.8,-0.5,-0.3, 0.4, 0.3,-0.2,-0.8},//ASP
        {-0.3,-0.4,-0.5,-0.7, 0.9,-0.7,-0.7,-0.4,-0.4,-0.3,-0.7,-0.7,-0.6,-0.6,-0.4, 0.0,-0.3,-0.8,-0.1,-0.3,-0.6,-0.7,-0.4,-0.8},//CYS
        {-0.1, 0.1, 0.0, 0.1,-0.7, 0.6, 0.2,-0.3, 0.3,-0.3,-0.2, 0.0,-0.1,-0.6, 0.0,-0.2,-0.2,-0.6,-0.5,-0.3, 0.0, 0.4,-0.1,-0.8},//GLN
        { 0.0,-0.3, 0.1, 0.3,-0.7, 0.2, 0.5,-0.1,-0.1,-0.3,-0.4,-0.1,-0.3,-0.7,-0.2,-0.1,-0.2,-0.8,-0.5,-0.3, 0.3, 0.4,-0.1,-0.8},//GLU
        { 0.1,-0.4, 0.0, 0.0,-0.4,-0.3,-0.1, 0.5,-0.4,-0.4,-0.5,-0.3,-0.4,-0.5,-0.2, 0.1,-0.1,-0.8,-0.6,-0.2, 0.0,-0.2,-0.2,-0.8},//GLY
        {-0.3, 0.1, 0.2, 0.0,-0.4, 0.3,-0.1,-0.4, 0.7,-0.4,-0.3,-0.2,-0.4,-0.3,-0.1,-0.2,-0.3,-0.3,-0.1,-0.3, 0.1, 0.1,-0.2,-0.8},//HIS
        {-0.1,-0.2,-0.2,-0.3,-0.3,-0.3,-0.3,-0.4,-0.4, 0.6, 0.1,-0.3, 0.1, 0.0,-0.3,-0.2, 0.0,-0.6,-0.2, 0.3,-0.3,-0.3,-0.1,-0.8},//ILE
        {-0.3,-0.4,-0.4,-0.5,-0.7,-0.2,-0.4,-0.5,-0.3, 0.1, 0.5,-0.4, 0.3, 0.0,-0.3,-0.4,-0.3,-0.3,-0.2, 0.1,-0.4,-0.3,-0.2,-0.8},//LEU
        {-0.2, 0.2, 0.1,-0.1,-0.7, 0.0,-0.1,-0.3,-0.2,-0.3,-0.4, 0.5, 0.0,-0.7,-0.2,-0.1,-0.1,-0.5,-0.5,-0.4, 0.0,-0.1,-0.2,-0.8},//LYS
        {-0.2,-0.1,-0.3,-0.4,-0.6,-0.1,-0.3,-0.4,-0.4, 0.1, 0.3, 0.0, 0.8,-0.1,-0.3,-0.2,-0.1,-0.6,-0.4, 0.1,-0.4,-0.2,-0.2,-0.8},//MET
        {-0.4,-0.5,-0.4,-0.7,-0.6,-0.6,-0.7,-0.5,-0.3, 0.0, 0.0,-0.7,-0.1, 0.8,-0.5,-0.3,-0.4,-0.1, 0.4,-0.3,-0.5,-0.6,-0.3,-0.8},//PHE
        { 0.1,-0.1,-0.2,-0.3,-0.4, 0.0,-0.2,-0.2,-0.1,-0.3,-0.3,-0.2,-0.3,-0.5, 0.6, 0.1,-0.1,-0.7,-0.6,-0.2,-0.2,-0.1,-0.2,-0.8},//PRO
        { 0.1,-0.1, 0.1, 0.0, 0.0,-0.2,-0.1, 0.1,-0.2,-0.2,-0.4,-0.1,-0.2,-0.3, 0.1, 0.3, 0.2,-0.2,-0.3,-0.2, 0.0,-0.1,-0.1,-0.8},//SER
        { 0.1,-0.2, 0.0,-0.1,-0.3,-0.2,-0.2,-0.1,-0.3, 0.0,-0.3,-0.1,-0.1,-0.4,-0.1, 0.2, 0.4,-0.6,-0.3, 0.0, 0.0,-0.2,-0.1,-0.8},//THR
        {-0.7, 0.1,-0.4,-0.8,-0.8,-0.6,-0.8,-0.8,-0.3,-0.6,-0.3,-0.5,-0.6,-0.1,-0.7,-0.2,-0.6, 1.2,-0.2,-0.8,-0.6,-0.7,-0.5,-0.8},//TRP
        {-0.4,-0.5,-0.2,-0.5,-0.1,-0.5,-0.5,-0.6,-0.1,-0.2,-0.2,-0.5,-0.4, 0.4,-0.6,-0.3,-0.3,-0.2, 0.8,-0.3,-0.3,-0.5,-0.3,-0.8},//TYR
        { 0.0,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3,-0.2,-0.3, 0.3, 0.1,-0.4, 0.1,-0.3,-0.2,-0.2, 0.0,-0.8,-0.3, 0.5,-0.3,-0.3,-0.1,-0.8},//VAL
        { 0.0,-0.2, 0.3, 0.4,-0.6, 0.0, 0.3, 0.0, 0.1,-0.3,-0.4, 0.0,-0.4,-0.5,-0.2, 0.0, 0.0,-0.6,-0.3,-0.3, 0.4, 0.2,-0.1,-0.8},//ASX
        {-0.1,-0.1, 0.0, 0.3,-0.7, 0.4, 0.4,-0.2, 0.1,-0.3,-0.3,-0.1,-0.2,-0.6,-0.1,-0.1,-0.2,-0.7,-0.5,-0.3, 0.2, 0.4,-0.1,-0.8},//GLX
        {-0.1,-0.2,-0.1,-0.2,-0.4,-0.1,-0.1,-0.2,-0.2,-0.1,-0.2,-0.2,-0.2,-0.3,-0.2,-0.1,-0.1,-0.5,-0.3,-0.1,-0.1,-0.1,-0.2,-0.8},//XXX
        {-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8, 0.1} //GAP
      },
      //PAM_160
      {
//        ALA  ARG  ASN  ASP  CYS  GLN  GLU  GLY  HIS  ILE  LEU  LYS  MET  PHE  PRO  SER  THR  TRP  TYR  VAL  ASX  GLX  XXX  GAP
        { 0.2,-0.2, 0.0, 0.0,-0.2,-0.1, 0.0, 0.1,-0.2,-0.1,-0.2,-0.2,-0.1,-0.3, 0.1, 0.1, 0.1,-0.5,-0.3, 0.0, 0.0, 0.0, 0.0,-0.7},//ALA
        {-0.2, 0.6,-0.1,-0.2,-0.3, 0.1,-0.2,-0.3, 0.1,-0.2,-0.3, 0.3,-0.1,-0.4,-0.1,-0.1,-0.1, 0.1,-0.4,-0.3,-0.1, 0.0,-0.1,-0.7},//ARG
        { 0.0,-0.1, 0.3, 0.2,-0.4, 0.0, 0.1, 0.0, 0.2,-0.2,-0.3, 0.1,-0.2,-0.3,-0.1, 0.1, 0.0,-0.4,-0.2,-0.2, 0.2, 0.1, 0.0,-0.7},//ASN
        { 0.0,-0.2, 0.2, 0.4,-0.5, 0.1, 0.3, 0.0, 0.0,-0.3,-0.4, 0.0,-0.3,-0.6,-0.2, 0.0,-0.1,-0.6,-0.4,-0.3, 0.3, 0.2,-0.1,-0.7},//ASP
        {-0.2,-0.3,-0.4,-0.5, 0.9,-0.5,-0.5,-0.3,-0.3,-0.2,-0.6,-0.5,-0.5,-0.5,-0.3, 0.0,-0.2,-0.7, 0.0,-0.2,-0.4,-0.5,-0.3,-0.7},//CYS
        {-0.1, 0.1, 0.0, 0.1,-0.5, 0.5, 0.2,-0.2, 0.2,-0.2,-0.2, 0.0,-0.1,-0.5, 0.0,-0.1,-0.1,-0.5,-0.4,-0.2, 0.1, 0.3,-0.1,-0.7},//GLN
        { 0.0,-0.2, 0.1, 0.3,-0.5, 0.2, 0.4, 0.0, 0.0,-0.2,-0.3,-0.1,-0.2,-0.5,-0.1, 0.0,-0.1,-0.7,-0.4,-0.2, 0.2, 0.3,-0.1,-0.7},//GLU
        { 0.1,-0.3, 0.0, 0.0,-0.3,-0.2, 0.0, 0.4,-0.3,-0.3,-0.4,-0.2,-0.3,-0.4,-0.1, 0.1,-0.1,-0.7,-0.5,-0.2, 0.0,-0.1,-0.1,-0.7},//GLY
        {-0.2, 0.1, 0.2, 0.0,-0.3, 0.2, 0.0,-0.3, 0.6,-0.3,-0.2,-0.1,-0.3,-0.2,-0.1,-0.1,-0.2,-0.3, 0.0,-0.2, 0.1, 0.1,-0.1,-0.7},//HIS
        {-0.1,-0.2,-0.2,-0.3,-0.2,-0.2,-0.2,-0.3,-0.3, 0.5, 0.2,-0.2, 0.2, 0.0,-0.2,-0.2, 0.0,-0.5,-0.2, 0.3,-0.2,-0.2,-0.1,-0.7},//ILE
        {-0.2,-0.3,-0.3,-0.4,-0.6,-0.2,-0.3,-0.4,-0.2, 0.2, 0.5,-0.3, 0.3, 0.1,-0.3,-0.3,-0.2,-0.2,-0.2, 0.1,-0.4,-0.3,-0.2,-0.7},//LEU
        {-0.2, 0.3, 0.1, 0.0,-0.5, 0.0,-0.1,-0.2,-0.1,-0.2,-0.3, 0.4, 0.0,-0.5,-0.2,-0.1, 0.0,-0.4,-0.4,-0.3, 0.0, 0.0,-0.1,-0.7},//LYS
        {-0.1,-0.1,-0.2,-0.3,-0.5,-0.1,-0.2,-0.3,-0.3, 0.2, 0.3, 0.0, 0.7, 0.0,-0.2,-0.2,-0.1,-0.4,-0.3, 0.1,-0.3,-0.2,-0.1,-0.7},//MET
        {-0.3,-0.4,-0.3,-0.6,-0.5,-0.5,-0.5,-0.4,-0.2, 0.0, 0.1,-0.5, 0.0, 0.7,-0.4,-0.3,-0.3,-0.1, 0.5,-0.2,-0.4,-0.5,-0.3,-0.7},//PHE
        { 0.1,-0.1,-0.1,-0.2,-0.3, 0.0,-0.1,-0.1,-0.1,-0.2,-0.3,-0.2,-0.2,-0.4, 0.5, 0.1, 0.0,-0.5,-0.5,-0.2,-0.1,-0.1,-0.1,-0.7},//PRO
        { 0.1,-0.1, 0.1, 0.0, 0.0,-0.1, 0.0, 0.1,-0.1,-0.2,-0.3,-0.1,-0.2,-0.3, 0.1, 0.2, 0.1,-0.2,-0.3,-0.1, 0.0,-0.1, 0.0,-0.7},//SER
        { 0.1,-0.1, 0.0,-0.1,-0.2,-0.1,-0.1,-0.1,-0.2, 0.0,-0.2, 0.0,-0.1,-0.3, 0.0, 0.1, 0.3,-0.5,-0.3, 0.0, 0.0,-0.1, 0.0,-0.7},//THR
        {-0.5, 0.1,-0.4,-0.6,-0.7,-0.5,-0.7,-0.7,-0.3,-0.5,-0.2,-0.4,-0.4,-0.1,-0.5,-0.2,-0.5, 1.2,-0.1,-0.6,-0.5,-0.6,-0.4,-0.7},//TRP
        {-0.3,-0.4,-0.2,-0.4, 0.0,-0.4,-0.4,-0.5, 0.0,-0.2,-0.2,-0.4,-0.3, 0.5,-0.5,-0.3,-0.3,-0.1, 0.8,-0.3,-0.3,-0.4,-0.3,-0.7},//TYR
        { 0.0,-0.3,-0.2,-0.3,-0.2,-0.2,-0.2,-0.2,-0.2, 0.3, 0.1,-0.3, 0.1,-0.2,-0.2,-0.1, 0.0,-0.6,-0.3, 0.4,-0.2,-0.2,-0.1,-0.7},//VAL
        { 0.0,-0.1, 0.2, 0.3,-0.4, 0.1, 0.2, 0.0, 0.1,-0.2,-0.4, 0.0,-0.3,-0.4,-0.1, 0.0, 0.0,-0.5,-0.3,-0.2, 0.3, 0.2,-0.1,-0.7},//ASX
        { 0.0, 0.0, 0.1, 0.2,-0.5, 0.3, 0.3,-0.1, 0.1,-0.2,-0.3, 0.0,-0.2,-0.5,-0.1,-0.1,-0.1,-0.6,-0.4,-0.2, 0.2, 0.3,-0.1,-0.7},//GLX
        { 0.0,-0.1, 0.0,-0.1,-0.3,-0.1,-0.1,-0.1,-0.1,-0.1,-0.2,-0.1,-0.1,-0.3,-0.1, 0.0, 0.0,-0.4,-0.3,-0.1,-0.1,-0.1,-0.1,-0.7},//XXX
        {-0.7,-0.7,-0.7,-0.7,-0.7,-0.7,-0.7,-0.7,-0.7,-0.7,-0.7,-0.7,-0.7,-0.7,-0.7,-0.7,-0.7,-0.7,-0.7,-0.7,-0.7,-0.7,-0.7, 0.1} //GAP
      },
      //PAM_250
      {
//        ALA  ARG  ASN  ASP  CYS  GLN  GLU  GLY  HIS  ILE  LEU  LYS  MET  PHE  PRO  SER  THR  TRP  TYR  VAL  ASX  GLX  XXX  GAP
        { 0.2,-0.2, 0.0, 0.0,-0.2, 0.0, 0.0, 0.1,-0.1,-0.1,-0.2,-0.1,-0.1,-0.3, 0.1, 0.1, 0.1,-0.6,-0.3, 0.0, 0.0, 0.0, 0.0,-0.8},//ALA
        {-0.2, 0.6, 0.0,-0.1,-0.4, 0.1,-0.1,-0.3, 0.2,-0.2,-0.3, 0.3, 0.0,-0.4, 0.0, 0.0,-0.1, 0.2,-0.4,-0.2,-0.1, 0.0,-0.1,-0.8},//ARG
        { 0.0, 0.0, 0.2, 0.2,-0.4, 0.1, 0.1, 0.0, 0.2,-0.2,-0.3, 0.1,-0.2,-0.3, 0.0, 0.1, 0.0,-0.4,-0.2,-0.2, 0.2, 0.1, 0.0,-0.8},//ASN
        { 0.0,-0.1, 0.2, 0.4,-0.5, 0.2, 0.3, 0.1, 0.1,-0.2,-0.4, 0.0,-0.3,-0.6,-0.1, 0.0, 0.0,-0.7,-0.4,-0.2, 0.3, 0.3,-0.1,-0.8},//ASP
        {-0.2,-0.4,-0.4,-0.5, 1.2,-0.5,-0.5,-0.3,-0.3,-0.2,-0.6,-0.5,-0.5,-0.4,-0.3, 0.0,-0.2,-0.8, 0.0,-0.2,-0.4,-0.5,-0.3,-0.8},//CYS
        { 0.0, 0.1, 0.1, 0.2,-0.5, 0.4, 0.2,-0.1, 0.3,-0.2,-0.2, 0.1,-0.1,-0.5, 0.0,-0.1,-0.1,-0.5,-0.4,-0.2, 0.1, 0.3,-0.1,-0.8},//GLN
        { 0.0,-0.1, 0.1, 0.3,-0.5, 0.2, 0.4, 0.0, 0.1,-0.2,-0.3, 0.0,-0.2,-0.5,-0.1, 0.0, 0.0,-0.7,-0.4,-0.2, 0.3, 0.3,-0.1,-0.8},//GLU
        { 0.1,-0.3, 0.0, 0.1,-0.3,-0.1, 0.0, 0.5,-0.2,-0.3,-0.4,-0.2,-0.3,-0.5, 0.0, 0.1, 0.0,-0.7,-0.5,-0.1, 0.0, 0.0,-0.1,-0.8},//GLY
        {-0.1, 0.2, 0.2, 0.1,-0.3, 0.3, 0.1,-0.2, 0.6,-0.2,-0.2, 0.0,-0.2,-0.2, 0.0,-0.1,-0.1,-0.3, 0.0,-0.2, 0.1, 0.2,-0.1,-0.8},//HIS
        {-0.1,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.3,-0.2, 0.5, 0.2,-0.2, 0.2, 0.1,-0.2,-0.1, 0.0,-0.5,-0.1, 0.4,-0.2,-0.2,-0.1,-0.8},//ILE
        {-0.2,-0.3,-0.3,-0.4,-0.6,-0.2,-0.3,-0.4,-0.2, 0.2, 0.6,-0.3, 0.4, 0.2,-0.3,-0.3,-0.2,-0.2,-0.1, 0.2,-0.3,-0.3,-0.1,-0.8},//LEU
        {-0.1, 0.3, 0.1, 0.0,-0.5, 0.1, 0.0,-0.2, 0.0,-0.2,-0.3, 0.5, 0.0,-0.5,-0.1, 0.0, 0.0,-0.3,-0.4,-0.2, 0.1, 0.0,-0.1,-0.8},//LYS
        {-0.1, 0.0,-0.2,-0.3,-0.5,-0.1,-0.2,-0.3,-0.2, 0.2, 0.4, 0.0, 0.6, 0.0,-0.2,-0.2,-0.1,-0.4,-0.2, 0.2,-0.2,-0.2,-0.1,-0.8},//MET
        {-0.3,-0.4,-0.3,-0.6,-0.4,-0.5,-0.5,-0.5,-0.2, 0.1, 0.2,-0.5, 0.0, 0.9,-0.5,-0.3,-0.3, 0.0, 0.7,-0.1,-0.4,-0.5,-0.2,-0.8},//PHE
        { 0.1, 0.0, 0.0,-0.1,-0.3, 0.0,-0.1, 0.0, 0.0,-0.2,-0.3,-0.1,-0.2,-0.5, 0.6, 0.1, 0.0,-0.6,-0.5,-0.1,-0.1, 0.0,-0.1,-0.8},//PRO
        { 0.1, 0.0, 0.1, 0.0, 0.0,-0.1, 0.0, 0.1,-0.1,-0.1,-0.3, 0.0,-0.2,-0.3, 0.1, 0.2, 0.1,-0.2,-0.3,-0.1, 0.0, 0.0, 0.0,-0.8},//SER
        { 0.1,-0.1, 0.0, 0.0,-0.2,-0.1, 0.0, 0.0,-0.1, 0.0,-0.2, 0.0,-0.1,-0.3, 0.0, 0.1, 0.3,-0.5,-0.3, 0.0, 0.0,-0.1, 0.0,-0.8},//THR
        {-0.6, 0.2,-0.4,-0.7,-0.8,-0.5,-0.7,-0.7,-0.3,-0.5,-0.2,-0.3,-0.4, 0.0,-0.6,-0.2,-0.5, 1.7, 0.0,-0.6,-0.5,-0.6,-0.4,-0.8},//TRP
        {-0.3,-0.4,-0.2,-0.4, 0.0,-0.4,-0.4,-0.5, 0.0,-0.1,-0.1,-0.4,-0.2, 0.7,-0.5,-0.3,-0.3, 0.0, 1.0,-0.2,-0.3,-0.4,-0.2,-0.8},//TYR
        { 0.0,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.1,-0.2, 0.4, 0.2,-0.2, 0.2,-0.1,-0.1,-0.1, 0.0,-0.6,-0.2, 0.4,-0.2,-0.2,-0.1,-0.8},//VAL
        { 0.0,-0.1, 0.2, 0.3,-0.4, 0.1, 0.3, 0.0, 0.1,-0.2,-0.3, 0.1,-0.2,-0.4,-0.1, 0.0, 0.0,-0.5,-0.3,-0.2, 0.3, 0.2,-0.1,-0.8},//ASX
        { 0.0, 0.0, 0.1, 0.3,-0.5, 0.3, 0.3, 0.0, 0.2,-0.2,-0.3, 0.0,-0.2,-0.5, 0.0, 0.0,-0.1,-0.6,-0.4,-0.2, 0.2, 0.3,-0.1,-0.8},//GLX
        { 0.0,-0.1, 0.0,-0.1,-0.3,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.2,-0.1, 0.0, 0.0,-0.4,-0.2,-0.1,-0.1,-0.1,-0.1,-0.8},//XXX
        {-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8, 0.1} //GAP
      }
    };

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a PAM_TABLE
    //! @param PAM_TABLE Pam table to be used
    AAAssignmentPAM::AAAssignmentPAM( const TableType PAM_TABLE) :
      m_PAMTable( PAM_TABLE)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to a new AAAssignmentPAM copied from this one
    AAAssignmentPAM *AAAssignmentPAM::Clone() const
    {
      return new AAAssignmentPAM( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AAAssignmentPAM::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the requested PAM table as a matrix
    //! @param PAM_TABLE requested PAM table name
    //! @return the requested PAM table as a matrix
    linal::Matrix< double> AAAssignmentPAM::GetPAMMatrix( const TableType PAM_TABLE)
    {
      // construct the matrix and return it
      return linal::Matrix< double>
      (
        biol::AATypes::s_NumberStandardAATypes + 4,
        biol::AATypes::s_NumberStandardAATypes + 4,
        &s_PAMTable[ PAM_TABLE][ 0][ 0]
      );
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that calculates the score between two assigned members
    //! @param MEMBER_A amino acid A that is compared
    //! @param MEMBER_B amino acid A that is compared
    //! @return value from the m_PAMTable for this combination of amino acids
    double AAAssignmentPAM::operator()( const biol::AABase &MEMBER_A, const biol::AABase &MEMBER_B) const
    {
      static const size_t matrix_dim( biol::AATypes::s_NumberStandardAATypes + 4);

      // check that the types are valid, if not, return undefined
      if( MEMBER_A.GetType() >= matrix_dim || MEMBER_B.GetType() >= matrix_dim)
      {
        return util::GetUndefined< double>();
      }

      // return matrix value for non-gaps
      return s_PAMTable[ m_PAMTable][ MEMBER_A.GetType()][ MEMBER_B.GetType()];
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write to ostream
    //! @param OSTREAM is the output stream
    //! @param INDENT indentation
    //! @return returns the output stream
    std::ostream &AAAssignmentPAM::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_PAMTable, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief read from istream
    //! @param ISTREAM is the input stream
    //! @return returns the input stream
    std::istream &AAAssignmentPAM::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_PAMTable, ISTREAM);

      // end
      return ISTREAM;
    }

  } // namespace score
} // namespace bcl

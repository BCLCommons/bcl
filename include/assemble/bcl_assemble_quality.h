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

#ifndef BCL_ASSEMBLE_QUALITY_H_
#define BCL_ASSEMBLE_QUALITY_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "align/bcl_align.fwd.hh"
#include "biol/bcl_biol.fwd.hh"
#include "io/bcl_io.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_compare.h"
#include "biol/bcl_biol_atom_types.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "quality/bcl_quality_measures.h"
#include "quality/bcl_quality_superimpose_measures.h"
#include "storage/bcl_storage_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Quality
    //! @brief convenience class for using different quality measure with Protein Models
    //! @details This class provides convenience functions for using different quality measures with Protein Models,
    //! in addition to calculation of RMSD100
    //! it is derived from FunctionInterface< ProteinModel, double> and is constructed with a quality measure, that is
    //! is employed in the operator, that calculates the quality measure for that protein to the native structure, which
    //! is stored in the aasequence in the chains
    //!
    //! @see @link example_assemble_quality.cpp @endlink
    //! @author woetzen
    //! @date Feb 3, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Quality :
      math::FunctionInterfaceSerializable< ProteinModel, double>
    {

    private:

    //////////
    // data //
    //////////

      //! the quality measure to be employed by the operator
      quality::Measure m_Measure;

      //! atom types to be considered for calculation
      storage::Set< biol::AtomType> m_AtomTypes;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief gives the flag which allows specifying multiple alignment files corresponding to multiple chains
      //! @return flag which allows specifying multiple alignment files corresponding to multiple chains
      static const util::ShPtr< command::FlagInterface> &GetFlagAlignments();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      //! @details construct with rmsd and backbone atom types
      Quality();

      //! @brief construct from a quality measure and the atom types to be considered
      //! @param QUALITY_MEASURE quality measure to be calculated
      //! @param ATOM_TYPES Set of atom types to be superimposed
      Quality
      (
        const quality::Measure &QUALITY_MEASURE,
        const storage::Set< biol::AtomType> &ATOM_TYPES
      );

      //! @brief virtual copy constructor
      //! @return pointer a new Quality copied from this quality
      Quality *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief calculate the given QUALITY_MEASURE between the coordinates in the SSEs of the Protein model and the original coordinates given from the pdbfile, that are stored in the original sequence
      //! @param QUALITY_MEASURE quality measure to be calculated
      //! @param PROTEIN_MODEL the protein model for which the QUALITY_MEASURE of the coordinates in the SSEs is calculated against the coordinates in the sequence in the chains
      //! @param ATOM_TYPES Set of atom types to be superimposed
      //! @return QUALITY_MEASURE between the coordinates in the SSEs of the Protein model and the original coordinates given from the pdbfile, that are stored in the original sequence
      static
      double Calculate
      (
        const quality::Measure &QUALITY_MEASURE,
        const ProteinModel &PROTEIN_MODEL,
        const storage::Set< biol::AtomType> &ATOM_TYPES
      );

      //! @brief calculate the given QUALITY_MEASURE between the coordinates in the SSEs of the Protein model and the original coordinates given from the pdbfile, that are stored in the original sequence
      //! @param QUALITY_MEASURE quality measure to be calculated
      //! @param PROTEIN_MODEL the protein model for which the QUALITY_MEASURE of the coordinates in the SSEs is calculated against the coordinates in the sequence in the chains
      //! @param NATIVE_MODEL native structure
      //! @param ATOM_TYPES Set of atom types to be superimposed
      //! @return QUALITY_MEASURE between the coordinates in the SSEs of the Protein model and the original coordinates given from the pdbfile, that are stored in the original sequence
      static
      double Calculate
      (
        const quality::Measure &QUALITY_MEASURE,
        const ProteinModel &PROTEIN_MODEL,
        const ProteinModel &NATIVE_MODEL,
        const storage::Set< biol::AtomType> &ATOM_TYPES
      );

      //! @brief superimpose the coordinates in the SSEs of the Protein model and the original coordinates given from the pdbfile, that are stored in the original sequence
      //! @param SUPERIMPOSE_MEASURE superimpose measure to be used for superimposition
      //! @param PROTEIN_MODEL protein model that will be superimpose to its native coordinates
      //! @param ATOM_TYPES Set of atom types to be superimposed
      //! @return whether superimposition succeeded and the transformation used
      static
      storage::Pair< bool, math::TransformationMatrix3D> SuperimposeModel
      (
        const quality::SuperimposeMeasure &SUPERIMPOSE_MEASURE,
        ProteinModel &PROTEIN_MODEL,
        const storage::Set< biol::AtomType> &ATOM_TYPES
      );

      //! @brief superimpose given model to another model
      //! @param SUPERIMPOSE_MEASURE superimpose measure to be used for superimposition
      //! @param PROTEIN_MODEL the protein model to be superimposed
      //! @param REFERENCE_MODEL reference structure
      //! @param ATOM_TYPES Set of atom types to be superimposed
      //! @return whether superimposition succeeded and the transformation used
      static
      storage::Pair< bool, math::TransformationMatrix3D> SuperimposeModel
      (
        const quality::SuperimposeMeasure &SUPERIMPOSE_MEASURE,
        ProteinModel &PROTEIN_MODEL,
        const ProteinModel &REFERENCE_MODEL,
        const storage::Set< biol::AtomType> &ATOM_TYPES
      );

      //! @brief superimpose given model to another model
      //! @param SUPERIMPOSE_MEASURE superimpose measure to be used for superimposition
      //! @param PROTEIN_MODEL the protein model to be superimposed
      //! @param REFERENCE_MODEL reference structure
      //! @param ATOM_TYPES Set of atom types to be superimposed
      //! @return whether superimposition succeeded and the transformation used
      static
      storage::Pair< bool, math::TransformationMatrix3D> SuperimposeModelWithAlignment
      (
        const quality::SuperimposeMeasure &SUPERIMPOSE_MEASURE,
        ProteinModel &PROTEIN_MODEL,
        const ProteinModel &REFERENCE_MODEL,
        const storage::Set< biol::AtomType> &ATOM_TYPES
      );

      //! @brief calculates RMSD100 by normalizing RMSD by the number of residues
      //! @param RMSD RMSD to be converted to RMSD100
      //! @param NUMBER_RESIDUES number of residues aligned in RMSD calculation
      //! @return RMSD100 calculated from the given RMSD
      static
      double
      RMSD100
      (
        const double RMSD,
        const size_t NUMBER_RESIDUES
      );

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that calculates quality measure of protein model to to coordinates stored in the aa sequences
      //! @param PROTEIN_MODEL the model to be considered
      //! @return the result of the quality measure this object was constructed with
      double operator()( const ProteinModel &PROTEIN_MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read ProteinModel from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write ProteinModel to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT nr of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    public:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief returns the pair of current and original coordinates for specified atomtypes for each amino acid that are in SSEs for the given chain
      //! @param CHAIN chain of interest
      //! @param ATOM_TYPES Set of atom types to be returned
      //! @return the pair of current and original coordinates for specified atomtypes for each amino acid that are in SSEs for the given chain
      static
      storage::Pair< util::SiPtrVector< const linal::Vector3D>, util::SiPtrVector< const linal::Vector3D> >
      GetSSEAndSequenceCoordinates
      (
        const Chain &CHAIN,
        const Chain &NATIVE,
        const storage::Set< biol::AtomType> &ATOM_TYPES
      );

      //! @brief returns the pair of current and original coordinates for specified atomtypes for each amino acid that are in SSEs for all the given chains
      //! @param MODEL Protein model of interest
      //! @param ATOM_TYPES Set of atom types to be returned
      //! @return the pair of current and original coordinates for specified atomtypes for each amino acid that are in SSEs for all the given chains
      static
      storage::Pair< util::SiPtrVector< const linal::Vector3D>, util::SiPtrVector< const linal::Vector3D> >
      GetChainAndSequenceCoordinates
      (
        const ProteinModel &MODEL,
        const storage::Set< biol::AtomType> &ATOM_TYPES
      );

      //! @brief GetAllDesiredAtomCoordinates puts all the desired atoms for both proteins into a StorageVectorND
      //! @param COMMON_AAS VectorND which holds the amino acids which are common to both proteins
      //! @param DESIRED_ATOMS Set which holds the atoms which are desired to be included in the calculation
      //! @return VectorND< 2, util::SiPtrVector< const linal::Vector3D> > which has desired atom coordinates
      static storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> > GetAllDesiredAtomCoordinates
      (
        const storage::VectorND< 2, util::SiPtrList< const biol::AABase> > &COMMON_AAS,
        const storage::Set< biol::AtomType> &DESIRED_ATOMS
      );

      //! @brief creates alignments from a single ProteinModel using the SSE residues
      //! @param PROTEIN_MODEL ProteinModel of interest
      //! @return alignments created from the ProteinModel using the SSE residues
      static storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > >
      CreateAlignmentFromProteinModelSSEs( const ProteinModel &PROTEIN_MODEL);

      //! @brief creates alignments from a single ProteinModel using specific residues identified by collector
      //! @param PROTEIN_MODEL ProteinModel of interest
      //! @param COLLECTOR the collector that has the specific residues that will be used
      //! @return alignments created from the ProteinModel using the specified residues
      static storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > >
      CreateAlignmentFromCollectorAASpecified
      (
        const ProteinModel &PROTEIN_MODEL, const CollectorAASpecified &COLLECTOR
      );

      //! @brief creates alignments from a single ProteinModel using the chain sequence residues
      //! @param PROTEIN_MODEL ProteinModel of interest
      //! @return alignments created from the ProteinModel using the chain sequence residues
      static storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > >
      CreateAlignmentFromProteinModelSequences( const ProteinModel &PROTEIN_MODEL);

      //! @brief creates alignments of chain sequences to SSE residues for a given ProteinModel
      //! @param PROTEIN_MODEL ProteinModel of interest
      //! @return alignments of chain sequences to SSE residues for a given ProteinModel
      static storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > >
      CreateAlignmentProteinModelSequencesToSSEs( const ProteinModel &PROTEIN_MODEL);

    public:

      //! @brief creates alignments of sequences in SSEs of two protein models
      //! @param PROTEIN_MODEL_A first ProteinModel of interest
      //! @param PROTEIN_MODEL_B second ProteinModel of interest
      //! @return alignments of sequences in SSEs of two protein models
      static storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > >
      CreateAlignmentProteinModels
      (
        const ProteinModel &PROTEIN_MODEL_A,
        const ProteinModel &PROTEIN_MODEL_B
      );

    protected:

      //! @brief aligns two single alignments (each containing 1 sequence) using the given comparison
      //! @param ALIGNMENT_A first alignment
      //! @param ALIGNMENT_B second alignment
      //! @param COMPARISON comparison object for two AABase
      //! @return combined alignments
      static storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > AlignSingleAlignments
      (
        storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > &ALIGNMENT_A,
        storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > &ALIGNMENT_B,
        const util::BinaryFunctionInterface< biol::AABase, biol::AABase, bool> &COMPARISON = biol::AACompareDataPtr()
      );

    public:

      //! @brief prunes the given alignment by removing all assignments with gaps
      //! @param ALIGNMENT Alignment of two protein sequences
      //! @return pruned alignment where the gaps were removed
      static align::AlignmentNode< biol::AABase> RemoveAssignmentsWithGapsFromAlignment
      (
        const align::AlignmentInterface< biol::AABase> &ALIGNMENT
      );

      //! @brief prunes the given alignment making sure each template amino acid has defined coordinates for given atom types
      //! @param ALIGNMENT Alignment of two protein sequences
      //! @param ATOM_TYPES Set of atom types for which the coordinates need to be defined
      //! @return pruned alignment where the residues where the template sequence has undefined coordinates were removed
      static align::AlignmentNode< biol::AABase> RemoveUndefinedTemplateAminoAcidsFromAlignment
      (
        const align::AlignmentInterface< biol::AABase> &ALIGNMENT,
        const storage::Set< biol::AtomType> &ATOM_TYPES
      );

      //! @brief prunes the given alignment making sure each amino acid has defined coordinates for given atom types
      //! @param ALIGNMENT Alignment of two protein sequences
      //! @param ATOM_TYPES Set of atom types of interest
      //! @return pruned alignment where the residues with undefined coordinates were removed
      static align::AlignmentNode< biol::AABase> RemoveUndefinedAminoAcidsFromAlignment
      (
        const align::AlignmentInterface< biol::AABase> &ALIGNMENT,
        const storage::Set< biol::AtomType> &ATOM_TYPES
      );

      //! @brief given an amino acid alignment and atom types generate two siptr vector of corresponding coordinates
      //!        for quality calculations
      //! @param ALIGNMENT Alignment of two protein sequences
      //! @param ATOM_TYPES Set of atom types of interest
      //! @return two siptr vector of corresponding coordinates
      static storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> > CoordinatesFromAlignment
      (
        const align::AlignmentInterface< biol::AABase> &ALIGNMENT,
        const storage::Set< biol::AtomType> &ATOM_TYPES
      );

      //! @brief given an amino acid alignment and atom types generate two siptr vector of corresponding coordinates
      //!        for quality calculations, and the number of residues involved (e.g. skipping undefined residues
      //! @param ALIGNMENT Alignment of two protein sequences
      //! @param ATOM_TYPES Set of atom types of interest
      //! @return two siptr vector of corresponding coordinates, and the number of residues in the alignment
      static storage::Pair< storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> >, size_t>
      CoordinatesAndResidueCountFromAlignment
      (
        const align::AlignmentInterface< biol::AABase> &ALIGNMENT,
        const storage::Set< biol::AtomType> &ATOM_TYPES
      );

      //! @brief given an amino acid alignments and atom types generate two siptr vector of corresponding coordinates
      //!        for quality calculations
      //! @param ALIGNMENTS Alignments of two protein sequences
      //! @param ATOM_TYPES Set of atom types of interest
      //! @return two siptr vector of corresponding coordinates
      static storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> > CoordinatesFromAlignments
      (
        const storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > &ALIGNMENTS,
        const storage::Set< biol::AtomType> &ATOM_TYPES
      );

      //! @brief uses a vector of filenames to read in alignments for chains - chain is last character before extension
      //! @param HANDLER the type of handler to read in the alignments
      //! @param ALIGN_FILES the files containing the alignments for each chain
      //! @return map with a character for the chain and a corresponding alignment for that chain
      static storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > GetChainAlignments
      (
        const align::HandlerInterface< biol::AABase> &HANDLER, const storage::Vector< io::DirectoryEntry> &ALIGN_FILES
      );

    }; // class Quality

  } // namespace assemble
} // namespace bcl

#endif //BCL_ASSEMBLE_QUALITY_H_

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

#ifndef BCL_ASSEMBLE_QUALITY_BATCH_H_
#define BCL_ASSEMBLE_QUALITY_BATCH_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "align/bcl_align.fwd.hh"
#include "biol/bcl_biol.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_atom_types.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "quality/bcl_quality_measures.h"
#include "storage/bcl_storage_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class QualityBatch
    //! @brief calculates qualities in a batch and returns a list of quality names and values
    //! @details This is function that iterates over member quality set and calculates the values and returns the
    //! results in a list of pairs of quality name and value. It also checks for any RMSD measures, and for these
    //! also calculates RMSD100 values and adds them to the result
    //!
    //! @see @link example_assemble_quality_batch.cpp @endlink
    //! @author karakam
    //! @date Nov 9, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API QualityBatch :
      public math::FunctionInterfaceSerializable< ProteinModel, storage::List< storage::Pair< std::string, double> > >
    {

    private:

    //////////
    // data //
    //////////

      //! set of qualities to calculate
      storage::Set< quality::Measure> m_Qualities;

      //! set of atoms to be used for quality calculations
      storage::Set< biol::AtomType> m_AtomTypes;

      //! suffix to add to end of quality names
      std::string m_QualitySuffix;

      //! whether to use readable names (for nice table output)
      bool m_UseReadableNames;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      QualityBatch();

      //! @brief constructor from set of qualities and set of atom types
      //! @param QUALITIES set of qualities to calculate
      //! @param ATOM_TYPES set of atom types to be used
      //! @param QUALITY_SUFFIX suffix to add to end of quality names
      //! @param USE_READABLE_NAMES use readable names (for nice table output)
      QualityBatch
      (
        const storage::Set< quality::Measure> &QUALITIES,
        const storage::Set< biol::AtomType> &ATOM_TYPES,
        const std::string &QUALITY_SUFFIX = "",
        const bool USE_READABLE_NAMES = false
      );

      //! @brief Clone function
      //! @return pointer to new QualityBatch
      QualityBatch *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief static function to create the names of columns to be generated from the qualities
      //! @param QUALITIES Set of qualities of interest
      //! @return vector of string that correspond to column names
      static storage::Vector< std::string> ColumnNamesFromQualities
      (
        const storage::Set< quality::Measure> &QUALITIES
      );

      //! @brief operator calculating the batch of qualities for the given model
      //! @param PROTEIN_MODEL protein model for which the qualities will be calculated
      //! @return list of pairs of quality names and values for the given model
      storage::List< storage::Pair< std::string, double> > operator()
      (
        const ProteinModel &PROTEIN_MODEL
      ) const;

      //! @brief operator calculating the batch of qualities for the two models
      //! @param PROTEIN_MODEL_A protein model for which the qualities will be calculated
      //! @param PROTEIN_MODEL_B protein model for which the qualities will be calculated
      //! @return list of pairs of quality names and values for the given models
      storage::List< storage::Pair< std::string, double> > operator()
      (
        const ProteinModel &PROTEIN_MODEL_A,
        const ProteinModel &PROTEIN_MODEL_B
      ) const;

      //! @brief construct a table that has protein distance measures to be append to a score table
      //! @param PROTEIN_MODEL to be used for measure calculations
      //! @return storage table of measures
      storage::Table< double> ConstructTable( const ProteinModel &PROTEIN_MODEL) const;

      //! @brief construct a table that has protein distance measures to be append to a score table
      //! @param PROTEIN_MODEL_A protein model for which the qualities will be calculated
      //! @param PROTEIN_MODEL_B protein model for which the qualities will be calculated
      //! @return storage table of measures
      storage::Table< double> ConstructTable
      (
        const ProteinModel &PROTEIN_MODEL_A,
        const ProteinModel &PROTEIN_MODEL_B
      ) const;

      //! @brief construct a table that has only the completeness estimate without checking for coordinates
      //! @param PROTEIN_MODEL to be used
      //! @return storage table with completeness estimate
      storage::Table< double> ContructTableWithCompletenessOnly( const ProteinModel &PROTEIN_MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief calculates a qualities from an alignment
      //! @param ALIGNMENTS the alignments from which the qualities will be calculated
      //! @return list of pairs of string indicating the quality and double indicating the quality value
      storage::List< storage::Pair< std::string, double> > QualitiesFromAlignment
      (
        const storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > &ALIGNMENTS,
        const storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > &SSE_ALIGNMENTS
      ) const;

      //! @brief construct a table that has protein distance measures to be append to a score table
      //! @param QUALITY_LIST list of calculated quality measures and their values
      //! @return storage table of measures
      storage::Table< double> ConstructTable
      (
        const storage::List< storage::Pair< std::string, double> > &QUALITY_LIST
      ) const;

      //! @brief get formatted name for better readability
      //! @param ROW_NAME original row name
      //! @return formatted name for better readability
      static std::string GetFormattedName( const std::string &ROW_NAME);

    }; // class QualityBatch

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_QUALITY_BATCH_H_

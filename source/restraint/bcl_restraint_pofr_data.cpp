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
#include "restraint/bcl_restraint_pofr_data.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_mutate_tree.h"
#include "fold/bcl_fold_score_weight_set.h"
#include "io/bcl_io_serialization.h"
#include "restraint/bcl_restraint_sas_debye.h"
#include "restraint/bcl_restraint_sas_density_data.h"
#include "restraint/bcl_restraint_sas_pofr.h"
#include "score/bcl_score_pofr.h"
#include "score/bcl_score_restraint_pofr.h"
#include "score/bcl_score_restraint_saxs.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    // initialize score
    fold::Score PofrData::e_ScorePofrRestraint( fold::GetScores().e_Undefined);

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> PofrData::s_Instance
    (
      util::Enumerated< Interface>::AddInstance( new PofrData())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from extension
    //! @param EXTENSION the extension used to identify files containing data for this restraint
    PofrData::PofrData() :
      m_Data(),
      m_Extension( GetDefaultExtension())
    {
    }

    //! @brief Clone function
    //! @return pointer to new PofrData
    PofrData *PofrData::Clone() const
    {
      return new PofrData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &PofrData::GetAlias() const
    {
      static const std::string s_name( "POFR");
      return s_name;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PofrData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the default file extension
    //! @return the default file extension
    const std::string &PofrData::GetDefaultExtension() const
    {
      static const std::string s_extension( ".pofr");
      return s_extension;
    }

    //! @brief gives reference to the specific type of data this restraint uses
    //! @return gives reference to the specific type of data this restraint uses
    const SasDensityData &PofrData::GetDensityData() const
    {
      return *m_Data;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer PofrData::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Interface derived class implements functionality for a restraint based on PofR data");
      serializer.AddInitializer
      (
        "Data",
        "experimental PofR curve",
        io::Serialization::GetAgent( &m_Data)
      );
      serializer.AddInitializer(
        "Extension",
        "the extension used to identify files containing pofr data",
        io::Serialization::GetAgent( &m_Extension)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief initialize the scores and add them to Scores enumerator
    void PofrData::InitializeScores()
    {
      if( !e_ScorePofrRestraint.IsDefined())
      {
        // read the restraints from the file, if they aren't already defined
        if( !m_Data.IsDefined())
        {
          // reset the data
          m_Data = util::Implementation< SasDensityData>( new SasDensityData());
          *m_Data = m_Data->ReadRestraintsFromFile();
        }

        SasPofR density_result;

        // Stores the experimental data in the Density Interface Class
        density_result.SetExperimentalDensity( m_Data);
        e_ScorePofrRestraint = fold::GetScores().AddScore
        (
          util::ShPtr< score::ProteinModel>
          (
            new score::RestraintPofr( density_result, score::PofR())
          )
        );
      }
    }

    //! @brief sets the weights of scores in a weight set
    //! @param SCORE_WEIGHT_SET the score weight set that will be modified
    void PofrData::ModifyScoreWeightSet( fold::ScoreWeightSet &SCORE_WEIGHT_SET) const
    {
      SCORE_WEIGHT_SET.SetWeight( e_ScorePofrRestraint, 5000);
    }

    //! @brief modify the mutate tree used
    //! @param MUTATE_TREE MutateTree to be modified
    void PofrData::ModifyMutateTree( fold::MutateTree &MUTATE_TREE) const
    {
      // set SSE remove to 0
      MUTATE_TREE.SetMutateTypeProbability( fold::MutateTree::e_Remove, 0);
    }

    //! @brief merges this protocol's mutate tree into given mutate tree
    //! @param MUTATE_TREE tree into which to merge this protocol's tree
    void PofrData::MergeAndModifyMutateTree( fold::MutateTree &MUTATE_TREE) const
    {
      util::ShPtr< fold::MutateTree> sp_mutate_tree( GetMutateTree());
      MUTATE_TREE.Merge( *sp_mutate_tree);
    }
 
  //////////////////////
  // input and output //
  //////////////////////

    //! @brief reads restraints formatted for this restraint type from an istream
    //! @return istream restraints formatted for this restraint type were read from
    std::istream &PofrData::ReadRestraints( std::istream &ISTREAM)
    {
      // reset the data
      m_Data = util::Implementation< SasDensityData>( new SasDensityData( m_Extension));

      // read from the stream
      m_Data->ReadFromDataFile( ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PofrData::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Data, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &PofrData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Data, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace restraint
} // namespace bcl

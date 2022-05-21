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
#include "restraint/bcl_restraint_saxs_data.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_mutate_tree.h"
#include "fold/bcl_fold_score_weight_set.h"
#include "restraint/bcl_restraint_sas_analysis.h"
#include "restraint/bcl_restraint_sas_debye.h"
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
    fold::Score SaxsData::e_ScoreSaxsRestraint( fold::GetScores().e_Undefined);

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SaxsData::s_Instance
    (
      util::Enumerated< Interface>::AddInstance( new SaxsData())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from extension
    //! @param EXTENSION the extension used to identify files containing data for this restraint
    SaxsData::SaxsData() :
      m_Data(),
      m_Extension( GetDefaultExtension())
    {
    }

    //! @brief Clone function
    //! @return pointer to new SaxsData
    SaxsData *SaxsData::Clone() const
    {
      return new SaxsData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &SaxsData::GetAlias() const
    {
      static const std::string s_name( "SAXS");
      return s_name;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SaxsData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the default file extension
    //! @return the default file extension
    const std::string &SaxsData::GetDefaultExtension() const
    {
      static const std::string s_extension( ".saxs");
      return s_extension;
    }

    //! @brief gives reference to the specific type of data this restraint uses
    //! @return gives reference to the specific type of data this restraint uses
    const SasScatteringData &SaxsData::GetScatteringData() const
    {
      return *m_Data;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief initialize the scores and add them to Scores enumerator
    void SaxsData::InitializeScores()
    {
      if( !e_ScoreSaxsRestraint.IsDefined())
      {
        // read the restraints from the file, if they aren't already defined
        if( !m_Data.IsDefined())
        {
          // reset the data
          m_Data = util::ShPtr< SasScatteringData>( new SasScatteringData());
          *m_Data = m_Data->ReadRestraintsFromFile();
        }

        bool approximate_loops( true);
        bool approximate_side_chains( true);
        double c1( 1.0);
        double c2( 0.0);
        bool cpu( false);
        bool sans( false);
        double deuterium_exchange( 0.0);

        // Setup Commandline Strings for either the opencl or non-opencl version of the code
        util::Implementation< SasDebyeInterface> saxs
        (
          SasAnalysis::SetDebyeImplementation
          (
            approximate_loops,
            approximate_side_chains,
            c1,
            c2,
            cpu,
            sans,
            deuterium_exchange
          )
        );

        //SasDebye saxs( true);
        //saxs.SetExperimentalData( m_Data);

        saxs->SetExperimentalData( m_Data);
        e_ScoreSaxsRestraint = fold::GetScores().AddScore
        (
          util::ShPtr< score::ProteinModel>
          (
            new score::RestraintSaxs( saxs, score::SasType())
          )
        );
      }
    }

    //! @brief sets the weights of scores in a weight set
    //! @param SCORE_WEIGHT_SET the score weight set that will be modified
    void SaxsData::ModifyScoreWeightSet( fold::ScoreWeightSet &SCORE_WEIGHT_SET) const
    {
      SCORE_WEIGHT_SET.SetWeight( e_ScoreSaxsRestraint, 500);
    }

    //! @brief modify the mutate tree used
    //! @param MUTATE_TREE MutateTree to be modified
    void SaxsData::ModifyMutateTree( fold::MutateTree &MUTATE_TREE) const
    {
      // set SSE remove to 0
      MUTATE_TREE.SetMutateTypeProbability( fold::MutateTree::e_Remove, 0);
    }

    //! @brief merges this protocol's mutate tree into given mutate tree
    //! @param MUTATE_TREE tree into which to merge this protocol's tree
    void SaxsData::MergeAndModifyMutateTree( fold::MutateTree &MUTATE_TREE) const
    {
      util::ShPtr< fold::MutateTree> sp_mutate_tree( GetMutateTree());
      MUTATE_TREE.Merge( *sp_mutate_tree);
    }
 
    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SaxsData::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription
      (
        "Allows use of small angle x-ray scattering data restraints. "
      );
      serializer.AddInitializer
      (
        "extension",
        "restraints will be read in from {path/prefix given by -restraint_prefix}{extension}",
        io::Serialization::GetAgent( &m_Extension),
        GetDefaultExtension()
      );
      return serializer;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief reads restraints formatted for this restraint type from an istream
    //! @return istream restraints formatted for this restraint type were read from
    std::istream &SaxsData::ReadRestraints( std::istream &ISTREAM)
    {
      // reset the data
      m_Data = util::ShPtr< SasScatteringData>( new SasScatteringData( m_Extension));

      // read from the stream
      m_Data->ReadFromDataFile( ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SaxsData::Read( std::istream &ISTREAM)
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
    std::ostream &SaxsData::Write( std::ostream &OSTREAM, const size_t INDENT) const
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

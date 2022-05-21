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

#ifndef BCL_APP_GENERATE_AA_PAIR_STATISTICS_H_
#define BCL_APP_GENERATE_AA_PAIR_STATISTICS_H_

// include the interface for all apps
#include "app/bcl_app_apps.h"

// include other forward headers - sorted alphabetically
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_interface.h"
#include "linal/bcl_linal_vector.h"
#include "math/bcl_math_running_average.h"
#include "sspred/bcl_sspred.h"
#include "storage/bcl_storage_vector.h"

namespace bcl
{
  namespace app
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class GenerateAAPairStatistics
    //! @brief Application for generating atom descriptors to determine the hybridization of different atom types and
    //!        determining covalent and van der waals radii from a molecule structure library, such as the cambridge
    //!        structural database
    //!
    //! @author mendenjl
    //! @date Aug 1, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API GenerateAAPairStatistics :
      public Interface
    {
    private:

    //////////
    // data //
    //////////

      //! Output filename base
      util::ShPtr< command::FlagInterface> m_OutputFilenameBase;

      //! Maximum distance to compute statistics for
      util::ShPtr< command::FlagInterface> m_MaxDistanceFlag;

      //! Initializer for protein model
      util::ShPtr< command::FlagInterface> m_InputFlag;

      //! SSPred type that should be used to determine the real statistics
      util::ShPtr< command::FlagInterface> m_SSPredFlag;

      //! Whether to calculate statistics for the membrane environment too
      util::ShPtr< command::FlagInterface> m_ConsiderEnvironmentFlag;

      //! Environment to calculate statistics for, if only one
      util::ShPtr< command::FlagInterface> m_RestrictEnvironmentFlag;

      mutable sspred::Method m_Method;                 //!< Method obtained from m_SSPredFlag
      mutable bool           m_ConsiderEnvironment;    //!< Whether to consider environment as well as secondary struct
      mutable size_t         m_MaxDistance;            //!< From m_MaxDistanceFlag
      mutable biol::EnvironmentType m_EnvironmentType; //!< Environment type of interest, if only considering one

      //! statistics, keyed by distance (outer vector), aa pair type [0-199] (inner vector)

      //! Counts of all times that the given AA pair are both in strands/helices
      mutable storage::Vector< linal::Vector< size_t> > m_StrandStrandCounts;
      mutable storage::Vector< linal::Vector< size_t> > m_HelixHelixCounts;
      mutable storage::Vector< linal::Vector< size_t> > m_CoilCoilCounts;
      mutable storage::Vector< linal::Vector< size_t> > m_MembraneMembraneCounts;
      mutable storage::Vector< linal::Vector< size_t> > m_TransitionTransitionCounts;
      mutable storage::Vector< linal::Vector< size_t> > m_SolubleSolubleCounts;

      //! Counts of all times that only one of the given AA pair is in a strand/helix
      mutable storage::Vector< linal::Vector< size_t> > m_StartsHelixCounts;
      mutable storage::Vector< linal::Vector< size_t> > m_StartsStrandCounts;
      mutable storage::Vector< linal::Vector< size_t> > m_StartsCoilCounts;
      mutable storage::Vector< linal::Vector< size_t> > m_EndsHelixCounts;
      mutable storage::Vector< linal::Vector< size_t> > m_EndsStrandCounts;
      mutable storage::Vector< linal::Vector< size_t> > m_EndsCoilCounts;
      mutable storage::Vector< linal::Vector< size_t> > m_StartsMembraneCounts;
      mutable storage::Vector< linal::Vector< size_t> > m_StartsTransitionCounts;
      mutable storage::Vector< linal::Vector< size_t> > m_StartsSolubleCounts;
      mutable storage::Vector< linal::Vector< size_t> > m_EndsMembraneCounts;
      mutable storage::Vector< linal::Vector< size_t> > m_EndsTransitionCounts;
      mutable storage::Vector< linal::Vector< size_t> > m_EndsSolubleCounts;

      //! Counts of occurrences of the AA pair type
      mutable storage::Vector< linal::Vector< size_t> > m_AAPairTypeCounts;

      //! Statistics for the euclidean distance between the two residues
      mutable storage::Vector< storage::Vector< math::RunningAverage< double> > > m_Separation;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      GenerateAAPairStatistics();

    public:

      // instantiate enumerator for GenerateAAPairStatistics class
      static const ApplicationType GenerateAAPairStatistics_Instance;

      //! @brief Clone function
      //! @return pointer to new GenerateAAPairStatistics
      GenerateAAPairStatistics *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief write the initializers for descriptor::AAPairStatistics
      //! @param OUTPUT output stream to write to
      void WriteStatisticsInitializer( std::ostream &OUTPUT) const;

      //! @brief compute the statistics requested
      //! @param MODEL the model of interest
      void Process( const assemble::ProteinModelWithCache &MODEL) const;

    }; // GenerateAAPairStatistics

  } // namespace app
} // namespace bcl

#endif // BCL_APP_GENERATE_AA_PAIR_STATISTICS_H_

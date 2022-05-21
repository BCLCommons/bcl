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
#include "biol/bcl_biol_ramachandran.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_types.h"
#include "biol/bcl_biol_environment_types.h"
#include "command/bcl_command_command_state.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "score/bcl_score.h"
#include "util/bcl_util_enumerated.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    //! @brief get the default static instance of this class
    //! @return static default instance of this class
    const Ramachandran &Ramachandran::GetDefaultInstance()
    {
      // initialize a static const instance
      static const Ramachandran s_default_instance( GetDefaultHistogramFilename());

      // end
      return s_default_instance;
    }

    //! @brief get the default SSTypeHistogramFilename
    //! @return the default SSTypeHistogramFilename
    const std::string &Ramachandran::GetDefaultHistogramFilename()
    {
      // initialize static const string to hold the default histogram filename
      static const std::string s_default_sstype_histogram_filename( "phi_psi_angles_by_sstype.histogram2D");

      // end
      return s_default_sstype_histogram_filename;
    }

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> Ramachandran::s_Instance
    (
      util::Enumerated< Ramachandran>::AddInstance( new Ramachandran())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Ramachandran::Ramachandran()
    {
    }

    //! @brief constructor from a AAType and SSType histogram filename
    //! @param SS_TYPE_HISTOGRAM_FILENAME filename for the phi/psi histogram according to SSTypes and AATypes
    Ramachandran::Ramachandran( const std::string &SS_TYPE_HISTOGRAM_FILENAME) :
      m_HistogramFilename( SS_TYPE_HISTOGRAM_FILENAME),
      m_AATypeMap(),
      m_SSTypeMap()
    {
      // initialize the members
      std::ostringstream oss;
      ReadInitializerSuccessHook( util::ObjectDataLabel(), oss);
    }

    //! @brief Clone function
    //! @return pointer to new Ramachandran
    Ramachandran *Ramachandran::Clone() const
    {
      return new Ramachandran( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Ramachandran::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the HistogramFilename
    //! @return the HistogramFilename
    const std::string &Ramachandran::GetHistogramFilename() const
    {
      return m_HistogramFilename;
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &Ramachandran::GetAlias() const
    {
      static const std::string s_alias( "Ramachandran");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer Ramachandran::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Distribution of phi/psi angles.");
      serializer.AddInitializer
      (
        "histogram filename",
        "path to the histogram file representing the distribution",
        io::Serialization::GetAgent( &m_HistogramFilename),
        GetDefaultHistogramFilename()
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief return a random phi psi for the given AAType
    //! @param AA_TYPE AAType of interest
    //! @return pair of random phi and psi values for the given AAType
    storage::VectorND< 2, double> Ramachandran::GetRandomPhiPsi( const AAType &AA_TYPE) const
    {
      // find the corresponding member in the AAType map
      storage::Map< AAType, HistogramDistributionPair>::const_iterator map_itr( m_AATypeMap.Find( AA_TYPE));

      // assert that it is found
      BCL_Assert( map_itr != m_AATypeMap.End(), "no distribution stored for AAType " + AA_TYPE.GetName());

      // get a random value
      return
        map_itr->second.Second().DetermineRandomCase2D
        (
          map_itr->second.First().GetBoundariesX().First(),
          map_itr->second.First().GetBoundariesY().First(),
          map_itr->second.First().GetBinSizeXY().First(),
          map_itr->second.First().GetBinSizeXY().Second()
        );
    }

    //! @brief return a random phi psi for the given AAType and SSType
    //! @param AA_TYPE AAType of interest
    //! @param SS_TYPE SSType of interest
    //! @return pair of random phi and psi values for the given AAType and SSType
    storage::VectorND< 2, double> Ramachandran::GetRandomPhiPsi
    (
      const AAType &AA_TYPE,
      const SSType SS_TYPE
    ) const
    {
      return this->GetRandomPhiPsi( AA_TYPE, SS_TYPE, GetEnvironmentTypes().e_Solution);
    }

    //! @brief return a random phi psi for the given AAType and SSType
    //! @param AA_TYPE AAType of interest
    //! @param SS_TYPE SSType of interest
    //! @return pair of random phi and psi values for the given AAType and SSType
    storage::VectorND< 2, double> Ramachandran::GetRandomPhiPsi
    (
      const AAType &AA_TYPE,
      const SSType SS_TYPE,
      const EnvironmentType &ENV_TYPE
    ) const
    {
      const storage::Map< SSType, storage::Map< AAType, HistogramDistributionPair> > &map
      (
        ENV_TYPE == GetEnvironmentTypes().e_MembraneCore ? m_SSTypeMapMembrane : m_SSTypeMap
      );

      // find the corresponding member in the AAType map
      storage::Map< SSType, storage::Map< AAType, HistogramDistributionPair> >::const_iterator
        map_itr_a( map.Find( SS_TYPE));

      // assert that it is found
      BCL_Assert( map_itr_a != map.End(), "no distribution stored for SSType " + SS_TYPE.GetName());

      // now find the distribution for the corresponding aatype
      storage::Map< AAType, HistogramDistributionPair>::const_iterator map_itr_b( map_itr_a->second.Find( AA_TYPE));

      // assert that a probability distribution is stored for this amino acid type
      BCL_Assert( map_itr_b != map_itr_a->second.End(), "no distribution stored for amino acid " + AA_TYPE.GetName());

      // get a random value
      return
        map_itr_b->second.Second().DetermineRandomCase2D
        (
          map_itr_b->second.First().GetBoundariesX().First(),
          map_itr_b->second.First().GetBoundariesY().First(),
          map_itr_b->second.First().GetBinSizeXY().First(),
          map_itr_b->second.First().GetBinSizeXY().Second()
        );
    }

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    //! @return return code indicating success or failure
    bool Ramachandran::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      if( !command::CommandState::IsInStaticInitialization())
      {
        Initialize();
      }
      return true;
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief function to read the ramachandran distributions and initialize the class members
    void Ramachandran::Initialize()
    {
      // reset
      m_SSTypeMap.Reset();
      m_AATypeMap.Reset();

      // initialize read
      io::IFStream read;

      // open read to the SSType PhiPsi histograms
      io::File::MustOpenIFStream( read, score::Score::AddHistogramPath( m_HistogramFilename));

      // iterate over the SSTypes
      for( size_t env_index( 0), n_env( 2); env_index < n_env; ++env_index)
      {
        std::string tmp;
        read >> tmp;
        for
        (
          SSTypes::const_iterator sstype_itr( GetSSTypes().Begin()), sstype_itr_end( GetSSTypes().COIL.GetIterator() + 1);
            sstype_itr != sstype_itr_end; ++sstype_itr
        )
        {
          // read in the enum
          SSType current_sstype;
          read >> current_sstype;

          // make sure it's correct
          BCL_Assert
          (
            *sstype_itr == current_sstype, "Wrong SSType " + sstype_itr->GetName() + " vs " + current_sstype.GetName()
          );

          // iterate over the valid AATypes
          for
          (
            AATypes::const_iterator
              aatype_itr( GetAATypes().Begin()),
              aatype_itr_end( GetAATypes().GetEnumIteratorFromIndex( AATypes::s_NumberStandardAATypes));
            aatype_itr != aatype_itr_end; ++aatype_itr
          )
          {
            // initialize pair to hold the histogram distribution pair
            HistogramDistributionPair hist_dist_pair;

            // read one letter code that represents the aatype
            std::string tmp;
            read >> tmp;

            // get the current amino acid type from the one letter code in "tmp"
            const AAType current_aatype( GetAATypes().AATypeFromOneLetterCode( tmp[ 0]));

            // assert that the aatypes in the histogram file are in the same order
            BCL_Assert
            (
              *aatype_itr == current_aatype,
              "unexpected aatype read from file! " + aatype_itr->GetName() + " != " +   current_aatype->GetName()
            );

            // read Histogram2D from stream
            read >> hist_dist_pair.First();

            if( !env_index)
            {
              // add to sstype independent map
              m_AATypeMap[ current_aatype].First().Combine( hist_dist_pair.First());
            }

            // create a Histogram2DDistribution from the histogram and store it
            hist_dist_pair.Second() = random::Histogram2DDistribution( hist_dist_pair.First());

            // insert into the member variable m_SSTypeMap
            ( env_index ? m_SSTypeMapMembrane : m_SSTypeMap)[ current_sstype][ current_aatype] = hist_dist_pair;
          }
        }
      }

      // clear the stream
      io::File::CloseClearFStream( read);

      // iterate over sstype independent map
      for
      (
        storage::Map< AAType, HistogramDistributionPair>::iterator
          itr( m_AATypeMap.Begin()), itr_end( m_AATypeMap.End());
        itr != itr_end;
        ++itr
      )
      {
        itr->second.Second() = random::Histogram2DDistribution( itr->second.First());
      }
    }

  } // namespace biol
} // namespace bcl

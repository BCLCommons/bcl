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
#include "descriptor/bcl_descriptor_aa_blast_profile.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "biol/bcl_biol_blast_profile_handler.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
  //////////
  // data //
  //////////

    // separate, anonymous namespace to prevent exporting symbols
    namespace
    {
      // add each of the possible instances to the enumerated instances
      util::ObjectInterface *AddInstances()
      {
        // keep a pointer to the last created instance
        util::ObjectInterface *last_instance( NULL);
        for( size_t method( 0); method < AABlastProfile::s_NumberMethods; ++method)
        {
          last_instance =
            util::Enumerated< Base< biol::AABase, float> >::AddInstance
            (
              new AABlastProfile( static_cast< AABlastProfile::Method>( method))
            );
        }
        return last_instance;
      }
    }

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AABlastProfile::s_Instance( AddInstances());

    //! @brief get the string for the method
    //! @param METHOD the method to retrieve the name for
    const std::string &AABlastProfile::GetMethodName( const Method &METHOD)
    {
      static const std::string s_names[ s_NumberMethods + 1] =
      {
        "AABlastProfile",
        "AABlastProbability",
        "AAType",
        "AA_BlastConservation",
        "AA_BlastPositionInformation",
        "AA_BlastAlignmentWeight",
        GetStaticClassName< Method>()
      };
      return s_names[ METHOD];
    }

    //! @brief get the description string for the method
    //! @param METHOD the method to retrieve the description for
    const std::string &AABlastProfile::GetMethodDescription( const Method &METHOD)
    {
      static const std::string s_names[ s_NumberMethods + 1] =
      {
        "BLAST-derived PSSM, 10 * Log10(Probability) for an AA",
        "BLAST-derived PSSM probability for an AA",
        "Binary descriptor (0/1), 1 iff the AA in question is the given type; order is ARNDCQEGHILKMFPSTWYV",
        "Reports the conservation of the amino acid = Sum(blast probability * log10 blast probability)",
        "Position information, as reported in the PSSM (second to last column)",
        "Alignment weight, relative to pseudocount weight, as reported in PSSM (last column)",
        GetStaticClassName< Method>()
      };
      return s_names[ METHOD];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor, allows choice of method
    AABlastProfile::AABlastProfile( const Method &METHOD) :
      m_Method( METHOD)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AABlastProfile
    AABlastProfile *AABlastProfile::Clone() const
    {
      return new AABlastProfile( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &AABlastProfile::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &AABlastProfile::GetAlias() const
    {
      return GetMethodName( m_Method);
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t AABlastProfile::GetNormalSizeOfFeatures() const
    {
      return int( m_Method) < int( s_FirstScalarMethod) ? biol::AATypes::s_NumberStandardAATypes : 1;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void AABlastProfile::Calculate
    (
      const iterate::Generic< const biol::AABase> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      if( m_Method == e_Identity)
      {
        STORAGE = float( 0.0);
        const size_t element_id( ELEMENT->GetType().GetIndex());
        if( element_id < biol::AATypes::s_NumberStandardAATypes)
        {
          STORAGE( element_id) = 1.0;
        }
        else if
        (
          ELEMENT->GetType().IsDefined()
          && ELEMENT->GetType()->GetParentType().GetIndex() < biol::AATypes::s_NumberStandardAATypes
        )
        {
          STORAGE( ELEMENT->GetType()->GetParentType().GetIndex()) = 1.0;
        }
        else
        {
          STORAGE = float( 1.0) / float( biol::AATypes::s_NumberStandardAATypes);
        }
        return;
      }

      if( !m_FileExtension.empty() && m_BlastProfiles.IsEmpty())
      {
        LoadFile();
      }

      if( m_BlastProfiles.IsEmpty() && !m_FileExtension.empty())
      {
        STORAGE = util::GetUndefined< float>();
        return;
      }

      // use the profile already stored on the aa if no file extension was given, otherwise, use the one loaded locally
      const biol::BlastProfile &profile
      (
        m_FileExtension.empty()
        ? ELEMENT->GetBlastProfile()
        : m_BlastProfiles( ELEMENT.GetPosition())
      );

      if( m_Method == e_Log10)
      {
        STORAGE.CopyValues( linal::Vector< float>( profile.GetProfile()));
      }
      else if( m_Method == e_Probability)
      {
        STORAGE.CopyValues( linal::Vector< float>( profile.GetProbabilities()));
      }
      else // scalar descriptors
      {
        switch( m_Method.GetEnum())
        {
          case e_Conservation: STORAGE( 0) = profile.CalculateConservation(); break;
          case e_Information: STORAGE( 0) = profile.GetPositionInformation(); break;
          case e_AlignmentWeight: STORAGE( 0) = profile.GetAlignmentWeight(); break;
          default:
            BCL_Exit( "Undefined enum!", -1);
        }
      }
    }

    //! @brief load the blast file
    void AABlastProfile::LoadFile()
    {
      if( m_FileExtension.empty() || m_Method == e_Identity)
      {
        return;
      }
      util::SiPtr< const assemble::ProteinModelWithCache> si_pmwc( this->GetCurrentObject());

      // Get the filename
      util::ShPtr< util::Wrapper< std::string> > sp_filename_wrapper
      (
        si_pmwc->GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile)
      );

      // Remove the last extension
      std::string basename
      (
        io::File::RemoveLastExtension( io::File::RemoveCompressionExtension( sp_filename_wrapper->GetData()))
      );

      m_BlastProfiles =
        biol::BlastProfileHandler::ReadProfilesForConstAASequence
        (
          *si_pmwc->GetChains().FirstElement()->GetSequence(),
          basename,
          m_FileExtension
        );

      if( m_BlastProfiles.GetSize() == si_pmwc->GetSize() || m_BlastProfiles.IsEmpty())
      {
        return;
      }

      // need to reorder the blast profiles, since not all the AAs are going to be iterated over
      storage::Vector< size_t> selected_aas;
      selected_aas.AllocateMemory( si_pmwc->GetSize());
      for( iterate::Generic< const biol::AABase> itr( si_pmwc->GetIterator()); itr.NotAtEnd(); ++itr)
      {
        selected_aas.PushBack( size_t( itr->GetSeqID() - 1));
      }
      m_BlastProfiles.Reorder( selected_aas);
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    //! In this case, resets the blast profiles
    void AABlastProfile::SetObjectHook()
    {
      m_BlastProfiles.Reset();
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AABlastProfile::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( GetMethodDescription( m_Method));
      if( m_Method != e_Identity)
      {
        serializer.AddOptionalInitializer
        (
          "extension",
          "ascii pssm file extension for the blast profile; if blank, uses the blast profile stored on the AA. "
          "When using dataset storage retrievers, this is specified by the blast extension parameter of "
          "ProteinDirectory or SequenceDirectory",
          io::Serialization::GetAgent( &m_FileExtension)
        );
      }
      return serializer;
    }

  } // namespace descriptor
} // namespace bcl

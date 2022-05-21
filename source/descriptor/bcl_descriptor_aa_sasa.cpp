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
#include "descriptor/bcl_descriptor_aa_sasa.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_count.h"
#include "assemble/bcl_assemble_aa_neighbor_list_container.h"
#include "assemble/bcl_assemble_aa_neighbor_list_container_generator_protein_model.h"
#include "assemble/bcl_assemble_aa_neighbor_vector.h"
#include "assemble/bcl_assemble_aa_sasa_ols.h"
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "util/bcl_util_enumerated.h"

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
        util::Enumerated< Base< biol::AABase, float> >::AddInstance
        (
          new AASasa
          (
            assemble::AASasaOLS(),
            "SASAOverlappingSpheres",
            "Computes SASA via the overlapping spheres method"
          )
        );
        util::Enumerated< Base< biol::AABase, float> >::AddInstance
        (
          new AASasa
          (
            assemble::AANeighborCount(),
            "SASANeighborCount",
            "Computes SASA using the count of AA neighbors"
          )
        );
        return
          util::Enumerated< Base< biol::AABase, float> >::AddInstance
          (
            new AASasa
            (
              assemble::AANeighborVector(),
              "SASANeighborVector",
              "Computes SASA using the neighbor vector method"
            )
          );
      }
    }

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AASasa::s_Instance( AddInstances());

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from exposure measure
    //! @param EXPOSURE_MEASURE exposure measure to use in order to calculate SASA
    //! @param ALIAS alias for this measure
    //! @param DESCRIPTION description of what this class does
    AASasa::AASasa
    (
      const assemble::AAExposureInterface &EXPOSURE_MEASURE,
      const std::string &ALIAS,
      const std::string &DESCRIPTION
    ) :
      m_AAExposure( EXPOSURE_MEASURE.Clone()),
      m_Alias( ALIAS),
      m_Description( DESCRIPTION + ", see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2712621/"),
      m_AANeighborListContainerGenerator
      (
        assemble::AANeighborListContainerGeneratorProteinModel::AANeighborListGenerator
        (
          m_AAExposure->GetDistanceCutoff(), m_AAExposure->GetMinimalSequenceSeparation(), true, false
        )
      )
    {
    }

    //! @brief Clone function
    //! @return pointer to new AASasa
    AASasa *AASasa::Clone() const
    {
      return new AASasa( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &AASasa::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &AASasa::GetAlias() const
    {
      return m_Alias;
    }

    //! @brief write neighbor counts and neighbor vectors
    //! @return no return
    void AASasa::WriteNCNV
    (
      const size_t MINIMAL_SEQUENCE_SEPARATION,
      std::ostream &OSTREAM,
      const char &CHAIN_ID,
      const assemble::ProteinModel &MODEL
    )
    {
      AASasa neighbor_count
      (
        assemble::AANeighborCount(),
        "SASANeighborCount",
        "Computes SASA using the count of AA neighbors"
      );

      // set minimal sequence separation to exclude residues very close in sequence
      neighbor_count.m_AAExposure->SetMinimalSequenceSeparation( MINIMAL_SEQUENCE_SEPARATION);
      neighbor_count.m_AANeighborListContainerGenerator =
        assemble::AANeighborListContainerGeneratorProteinModel::AANeighborListGenerator
        (
          neighbor_count.m_AAExposure->GetDistanceCutoff(), neighbor_count.m_AAExposure->GetMinimalSequenceSeparation(), true, false
        );

      AASasa neighbor_vector
      (
        assemble::AANeighborVector(),
        "SASANeighborVector",
        "Computes SASA using the neighbor vector method"
      );

      // create a protein model with cache so that it can be used to calculate descriptors
      assemble::ProteinModelWithCache pmwc( MODEL, true);
      neighbor_count.SetObject( pmwc);
      neighbor_vector.SetObject( pmwc);

      OSTREAM << "# SASA Neighbor Counts and Neighbor Vector Calculation via BCL Algorithm on a structure with " << MODEL.GetNumberOfChains() << " chains\n";
      OSTREAM << "# SeqId\tNeighborCount\tNeighborVector\n";
      for( Iterator< biol::AABase> itr( pmwc.GetIterator()); itr.NotAtEnd(); ++itr)
      {
        // skip the aa if its not the right chain id
        if( itr( 0)->GetChainID() != CHAIN_ID)
        {
          continue;
        }

        const biol::AABase &aa( *itr( 0));

        // write out the aa identification, followed by sasa nc and nv
        OSTREAM << util::Format().W( 5).R()( aa.GetSeqID()) << ',' << neighbor_count( itr)( 0) << ',' << neighbor_vector( itr)( 0) << '\n';
      }
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t AASasa::GetNormalSizeOfFeatures() const
    {
      return size_t( 1);
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void AASasa::Calculate
    (
      const iterate::Generic< const biol::AABase> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      // find the entry for given amino acid
      const storage::Map< util::SiPtr< const biol::AABase>, float, biol::AALessThanSeqID>::const_iterator itr
      (
        m_AASasaStorage.Find( *ELEMENT)
      );

      // check if the sasa for this aa exists
      if( itr == m_AASasaStorage.End())
      {
        // no sasa for this AA exists, set storage to undefined
        STORAGE = util::GetUndefined< float>();
      }
      else
      {
        // no sasa for this AA exists, return undefined
        STORAGE( 0) = itr->second;
      }
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AASasa::GetSerializer() const
    {
      return io::Serializer().SetClassDescription( m_Description);
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void AASasa::SetObjectHook()
    {
      // reset the map to store the AA's for the new protein model
      m_AASasaStorage.Reset();

      // get a pointer to the protein
      util::SiPtr< const assemble::ProteinModel> sp_protein_model( this->GetCurrentObject());

      // calculate SASA on the whole protein model

      // collect AANeighborLists for all amino acids in the given protein model
      const assemble::AANeighborListContainer all_aa_neighbor_list
      (
        m_AANeighborListContainerGenerator->operator ()( *sp_protein_model)
      );

      // iterate over all lists in the all_aa_neighbor_list
      for
      (
        assemble::AANeighborListContainer::const_iterator
          aa_itr( all_aa_neighbor_list.Begin()), aa_itr_end( all_aa_neighbor_list.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        m_AASasaStorage[ aa_itr->first] = m_AAExposure->operator ()( aa_itr->second);
      }
    }

  } // namespace descriptor
} // namespace bcl

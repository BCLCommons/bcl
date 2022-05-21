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
#include "descriptor/bcl_descriptor_aa_phi_psi.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_serialization.h"
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
        util::Enumerated< Base< biol::AABase, float> >::AddInstance( new AAPhiPsi( true, -10.0));
        return
          util::Enumerated< Base< biol::AABase, float> >::AddInstance( new AAPhiPsi( false, -90.0));
      }
    }

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AAPhiPsi::s_Instance( AddInstances());

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from angle type and origin
    //! @param PHI true to calculate phi angles (else, calculate psi)
    //! @param ORIGIN origin point; wrap angles below this around by adding 2 pi to them
    AAPhiPsi::AAPhiPsi( const bool &PHI, const float &ORIGIN) :
      m_CalculatePhi( PHI),
      m_Origin( ORIGIN)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AAPhiPsi
    AAPhiPsi *AAPhiPsi::Clone() const
    {
      return new AAPhiPsi( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &AAPhiPsi::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &AAPhiPsi::GetAlias() const
    {
      static const std::string s_phi( "Phi"), s_psi( "Psi");
      return m_CalculatePhi ? s_phi : s_psi;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t AAPhiPsi::GetNormalSizeOfFeatures() const
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
    void AAPhiPsi::Calculate
    (
      const iterate::Generic< const biol::AABase> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      // find the entry for given amino acid
      const storage::Map< util::SiPtr< const biol::AABase>, float, biol::AALessThanSeqID>::const_iterator itr
      (
        m_AAPhiPsiStorage.Find( *ELEMENT)
      );

      // check if the sasa for this aa exists
      if( itr == m_AAPhiPsiStorage.End() || !itr->first.IsDefined())
      {
        // no phi/psi for this AA exists, set storage to undefined
        STORAGE( 0) = util::GetUndefined< float>();
      }
      else
      {
        // return the calculated phi or psi angle
        STORAGE( 0) = itr->second;
      }
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AAPhiPsi::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        m_CalculatePhi
        ? "Calculates the phi angle (C-C-N-C)"
        : "Calculates the psi angle (N-C-C-N)"
      );
      parameters.AddInitializer
      (
        "origin",
        "point at which to wrap all other angles around by adding 2-pi",
        io::Serialization::GetAgentWithRange( &m_Origin, -360.0, 0.0),
        util::Format()( m_Origin)
      );
      return parameters;
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void AAPhiPsi::SetObjectHook()
    {
      // reset the map to store the AA's for the new protein model
      m_AAPhiPsiStorage.Reset();

      // get a pointer to the protein
      util::SiPtr< const assemble::ProteinModel> sp_protein_model( this->GetCurrentObject());

      // calculate phi-psi angles for the whole protein model
      util::SiPtrVector< const biol::AASequence> sequences( sp_protein_model->GetSequences());
      for
      (
        util::SiPtrVector< const biol::AASequence>::const_iterator
          itr_seq( sequences.Begin()), itr_seq_end( sequences.End());
        itr_seq != itr_seq_end;
        ++itr_seq
      )
      {
        const biol::AASequence &seq( **itr_seq);
        if( seq.GetSize() <= size_t( 1))
        {
          continue;
        }
        if( m_CalculatePhi)
        {
          // phi
          for
          (
            biol::AASequence::const_iterator
              itr_prev( seq.Begin()), itr( itr_prev + 1), itr_end( seq.End());
            itr != itr_end;
            ++itr, ++itr_prev
          )
          {
            // retrieve the carbon of previous amino acid
            const biol::Atom &previous_c( ( *itr_prev)->GetAtom( biol::GetAtomTypes().C));
            float phi( math::Angle::Degree( ( *itr)->CalculatePhi( previous_c)));
            // handle wrap-around
            if( phi < m_Origin)
            {
              phi += 360.0;
            }
            else if( phi >= m_Origin + 360.0)
            {
              phi -= 360.0;
            }
            m_AAPhiPsiStorage[ **itr] = phi;
          }
        }
        else
        {
          // psi
          for
          (
            biol::AASequence::const_iterator
              itr_prev( seq.Begin()), itr( itr_prev + 1), itr_end( seq.End());
            itr != itr_end;
            ++itr, ++itr_prev
          )
          {
            // retrieve the nitrogen of the next amino acid
            const biol::Atom &next_n( ( *itr)->GetAtom( biol::GetAtomTypes().N));
            float psi( math::Angle::Degree( ( *itr_prev)->CalculatePsi( next_n)));
            // handle wrap-around
            if( psi < m_Origin)
            {
              psi += 360.0;
            }
            else if( psi >= m_Origin + 360.0)
            {
              psi -= 360.0;
            }
            m_AAPhiPsiStorage[ **itr_prev] = psi;
          }
        }
      }
    }

  } // namespace descriptor
} // namespace bcl

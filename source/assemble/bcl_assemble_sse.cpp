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
#include "assemble/bcl_assemble_sse.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse_compare.h"
#include "assemble/bcl_assemble_sse_geometry.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "biol/bcl_biol_aa_sequence_phi_psi.h"
#include "quality/bcl_quality_rmsd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> SSE::s_Instance
    (
      GetObjectInstances().AddInstance( new SSE())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from SSType
    //! @param SS_TYPE specific SSType
    SSE::SSE( const biol::SSType &SS_TYPE) :
      biol::AASequence(),
      m_SSType( SS_TYPE),
      m_Orientation(),
      m_XExtent( m_SSType->GetRadialExtent()),
      m_YExtent( m_SSType->GetRadialExtent()),
      m_Fragments()
    {
    }

    //! @brief construct by AASequence and SSType
    //! @param SEQUENCE amino acid sequence
    //! @param SS_TYPE specific SSType
    SSE::SSE( const biol::AASequence &SEQUENCE, const biol::SSType &SS_TYPE) :
      biol::AASequence( SEQUENCE),
      m_SSType( SS_TYPE),
      m_Orientation(),
      m_XExtent( m_SSType->GetRadialExtent()),
      m_YExtent( m_SSType->GetRadialExtent()),
      m_Fragments()
    {
      // only set the geometry if the sequence has proper coordinates
      if( SEQUENCE.HasDefinedCoordinates())
      {
        SetGeometry();
      }
    }

    //! @brief copy constructor
    //! @param SSE_TO_COPY
    SSE::SSE( const SSE &SSE_TO_COPY) :
      biol::AASequence( SSE_TO_COPY),
      m_SSType( SSE_TO_COPY.m_SSType),
      m_Orientation( SSE_TO_COPY.m_Orientation),
      m_XExtent( SSE_TO_COPY.m_XExtent),
      m_YExtent( SSE_TO_COPY.m_YExtent),
      m_Fragments( SSE_TO_COPY.m_Fragments.HardCopy())
    {
    }

    //! @brief Clone function
    //! @return pointer to new SSE
    SSE *SSE::Clone() const
    {
      return new SSE( *this);
    }

    //! @brief virtual hard copy all amino acids and their data
    //! @return SSE with independent hard copied AADatas
    SSE *SSE::HardCopy() const
    {
      // copy for this sequence
      SSE *new_sse( Clone());

      // iterate over all amino acids to make hard copy of the Data
      for
      (
        biol::AASequence::iterator aa_itr( new_sse->Begin()), aa_itr_end( new_sse->End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // hard copy the data for that aa
        ( *aa_itr)->SetData( ( *aa_itr)->GetData().HardCopy());
      }

      // end
      return new_sse;
    }

    //! @brief destructor
    SSE::~SSE()
    {
      m_GeometryDestructorSignal.Emit( *this);
      m_DestructorSignal.Emit( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SSE::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get identification of this SSE
    //! @return string with identification for this SSE
    std::string SSE::GetIdentification() const
    {
      // initialize identification with SSType, seqid and 3 letter code for beginning and ending amino acid
      std::string identification
      (
        m_SSType.GetName() + ' ' + util::Format()( GetChainID())
      );

      // add identifiers for the first and last residues this sse is spanning
      if( GetSize() > 0)
      {
        identification += ' ' +
          util::Format().W( 4)( GetFirstAA()->GetSeqID()) + ' ' + GetFirstAA()->GetType()->GetThreeLetterCode() +
          " <==> " +
          util::Format().W( 4)( GetLastAA()->GetSeqID()) + ' ' + GetLastAA()->GetType()->GetThreeLetterCode();
      }

      // end
      return identification;
    }

    //! @brief returns the requested extent
    //! @param AXIS axis of interest
    //! @return the requested extent
    double SSE::GetExtent( const coord::Axis &AXIS) const
    {
      // return the appropriate extent
      if( AXIS == coord::GetAxes().e_X)
      {
        return m_XExtent;
      }
      if( AXIS == coord::GetAxes().e_Y)
      {
        return m_YExtent;
      }
      if( AXIS == coord::GetAxes().e_Z)
      {
        return GetLength() / 2.0;
      }

      // else return undefined
      return util::GetUndefinedDouble();
    }

    //! @brief sets the extents
    //! @param EXTENT vector containing the extents in x, y, z
    void SSE::SetExtents( const linal::Vector3D &EXTENT)
    {
      m_XExtent = EXTENT.X();
      m_YExtent = EXTENT.Y();
    }

    //! @brief sets the SSType
    //! @param SS_TYPE SSType to be set
    void SSE::SetType( const biol::SSType &SS_TYPE)
    {
      m_SSType = SS_TYPE;
      SetExtents();
    }

    //! @brief Get SSE hash string to aid in identifying similar chains
    std::string SSE::GetHashString() const
    {
      return util::Format()( m_SSType->GetOneLetterCode())
             + ':' + util::Format()( this->GetFirstAA()->GetSeqID())
             + '-' + util::Format()( this->GetLastAA()->GetSeqID());
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief translate the object along a given TRANSLATION vector
    //! @param TRANSLATION Translation to be applied
    void SSE::Translate( const linal::Vector3D &TRANSLATION)
    {
      // translate sequence
      AASequence::Translate( TRANSLATION);

      // create a transformation matrix from the translation
      math::TransformationMatrix3D transform;
      transform( TRANSLATION);

      // translate the geometries
      TransformGeometries( transform);
    }

    //! @brief transform the object by a given TransformationMatrix3D
    //! @param TRANSFORMATION_MATRIX_3D TransformationMatrix3D to be applied
    void SSE::Transform( const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D)
    {
      // transform sequence
      AASequence::Transform( TRANSFORMATION_MATRIX_3D);

      // transform geometries
      TransformGeometries( TRANSFORMATION_MATRIX_3D);
    }

    //! @brief rotate the object by a given RotationMatrix3D
    //! @param ROTATION_MATRIX_3D RotationMatrix3D to be applied
    void SSE::Rotate( const math::RotationMatrix3D &ROTATION_MATRIX_3D)
    {
      // rotate sequence
      AASequence::Rotate( ROTATION_MATRIX_3D);

      // rotate geometries
      TransformGeometries( math::TransformationMatrix3D( ROTATION_MATRIX_3D));
    }

    //! @brief rotation around ROTATION_AXIS that goes through ROTATION_POINT and given ANGLE
    //! @param ROTATION_POINT point which rotation should pass through
    //! @param ROTATION_AXIS RotationAxis around which SSE will be rotated
    //! @param ANGLE angle for the rotation
    void SSE::Rotate( const linal::Vector3D &ROTATION_POINT, const linal::Vector3D &ROTATION_AXIS, const double &ANGLE)
    {
      if( ANGLE == 0 || ( ROTATION_AXIS.X() == 0.0 && ROTATION_AXIS.Y() == 0.0 && ROTATION_AXIS.Z() == 0.0))
      {
        BCL_MessageCrt( "no rotation");

        return;
      }

      // translate all coordinates to origin
      math::TransformationMatrix3D matrix( -ROTATION_POINT);

      //build TransformationMatrix3D
      matrix( math::RotationMatrix3D( ROTATION_AXIS, ANGLE));

      // translate all coordinates to former position
      matrix( ROTATION_POINT);

      // transform with TransformationMatrix
      Transform( matrix);
    }

    //! @brief sets the conformation to idealized at origin
    void SSE::SetToIdealConformationAtOrigin()
    {
      if( !m_SSType->IsStructured())
      {
        return;
      }

      // idealize this sequence
      biol::AASequenceFactory::IdealizeSequence( *this, m_SSType);

      // reset the orientation
      m_Orientation = math::TransformationMatrix3D();

      // iterate over the geometries
      m_Fragments.Reset();
      SetFragmentGeometries();
    }

    //! @brief sets the conformation to idealized in place
    void SSE::SetToIdealConformationInPlace()
    {
      if( !m_SSType->IsStructured())
      {
        return;
      }

      // make a copy of the sequence
      biol::AASequence template_copy( *this);

      // reset the fragments
      m_Fragments.Reset();

      // set sequence to ideal conformation
      SetToIdealConformationAtOrigin();

      // superimpose with template
      Transform( biol::AASequenceFactory::CalculateSuperimposition( *this, template_copy));

      // set the geometries
      SetFragmentGeometries();
    }

    //! @brief sets all geometries; main geometry and fragment geometries
    void SSE::SetGeometry()
    {
      // set the main geometry
      SetMainGeometry();

      // set fragment geometries
      SetFragmentGeometries();
    }

    //! @brief sets the main geometry to the geometry of an copy SSE that is idealized, does not set fragment geometries
    void SSE::SetMainGeometry()
    {
      // set the extents
      SetExtents();

      // if this is a coil
      if( !m_SSType->IsStructured())
      {
        // initialize a transformation matrix
        // calculate the center of mass
        m_Orientation =
        math::TransformationMatrix3D
        (
          *coord::GetAxes().e_X,
          *coord::GetAxes().e_Y,
          *coord::GetAxes().e_Z,
          coord::CenterOfMass( GetAtomCoordinates( biol::GetAtomTypes().GetBackBoneAtomTypes()), true)
        );
      }

      // else if it's a helix or a strand
      else
      {
        // make a copy of the sequence
        biol::AASequence sequence_copy( *this);

        // idealize the copy
        biol::AASequenceFactory::IdealizeSequence( sequence_copy, m_SSType);

        // calculate the transformation that will superimpose the idealized sequence onto this sequence
        math::TransformationMatrix3D transformation
        (
          biol::AASequenceFactory::CalculateSuperimposition( sequence_copy, *this)
        );

        // update the main geometry to this transformation
        m_Orientation = transformation;
      }
    }

    //! @brief returns a storage::Set of SSEs chopped to pieces of minimal size with spacing of one aa
    //! @param SIZE minimal size
    //! @return a storage::Set of SSEs chopped to pieces of minimal size with spacing of one aa
    storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap> SSE::Chop( const size_t &SIZE) const
    {
      BCL_Assert( SIZE != 0, "cannot chop into sequences of size 0");

      //chop the sequence
      util::ShPtrVector< biol::AASequence> chopped_sequence( biol::AASequence::ChopSequence( SIZE));

      //build new sses from chopped sequence
      storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap> chopped_sses;

      //iterate over chopped sequence and pushback new sses
      for
      (
        util::ShPtrVector< biol::AASequence>::const_iterator seq_itr( chopped_sequence.Begin()),
          seq_itr_end( chopped_sequence.End());
        seq_itr != seq_itr_end;
        ++seq_itr
      )
      {
        chopped_sses.Insert( util::ShPtr< SSE>( new SSE( **seq_itr, GetType())));
      }

      //return chopped sequence
      return chopped_sses;
    }

    //! @brief finds the coordinates for the given amino acid and prepends it to the sequence
    //! @param AMINO_ACID amino acid to be prepended
    //! @param IDEALIZE bool whether to idealize the SSE prior to prepending
    void SSE::Prepend( const biol::AABase &AMINO_ACID, const bool IDEALIZE)
    {
      // check the chain id and the sequence id of the given amino acid for compatibility
      BCL_Assert
      (
        AMINO_ACID.GetChainID() == GetChainID() && AMINO_ACID.GetSeqID() == GetFirstAA()->GetSeqID() - 1,
        "The given amino acid: " + AMINO_ACID.GetIdentification() + " cannot be added to this SSE " +
        GetIdentification()
      );

      // initialize transformation matrix
      math::TransformationMatrix3D transform;

      // if the SSE should be idealized
      if( IDEALIZE)
      {
        // idealize this SSE
        SetToIdealConformationInPlace();

        // initialize the transformation matrix to the body of the fragment
        transform = this->GetOrientation();

        // find the translation
        transform( -1.0 * GetLength() / 2.0 * GetAxis( coord::GetAxes().e_Z));
      }
      // do not idealize SSE and determine transformation from last residue
      else
      {
        // get the first AA
        util::ShPtr< biol::AABase> first_aa( GetData().FirstElement());

        // make a hardcopy and set to ideal at origin
        util::ShPtr< biol::AABase> copy_aa( first_aa.HardCopy());
        copy_aa->SetToIdealConformation( m_SSType, math::TransformationMatrix3D());

        // get the transformation matrix to superimpose the aa's
        transform = quality::RMSD::SuperimposeCoordinates( first_aa->GetAtomCoordinates(), copy_aa->GetAtomCoordinates());
      }

      // make a copy for the transformation for going from one residue to the next and invert it
      math::TransformationMatrix3D transform_for_residue( m_SSType->GetTransformationMatrixForResidues());
      transform_for_residue.Invert();

      // make a copy of this amino acid
      util::ShPtr< biol::AABase> new_aa( AMINO_ACID.Clone());

      // now set to ideal coordinates for the calculated transformation at the origin
      new_aa->SetToIdealConformation( m_SSType, transform_for_residue);

      // now to transform to the original coordinate
      new_aa->Transform( transform);

      // prepend the amino acid into the sequence
      biol::AASequence::PushFront( *new_aa);

      // set geometry
      SetGeometry();
    }

    //! @brief calculates ideal coordinates for the given amino acid and appends it to the sequence
    //! @param AMINO_ACID amino acid to be appended
    //! @param IDEALIZE bool whether to idealize the SSE prior to appending
    void SSE::Append( const biol::AABase &AMINO_ACID, const bool IDEALIZE)
    {
      // check the chain id and the sequence id of the given amino acid for compatibility
      BCL_Assert
      (
        AMINO_ACID.GetChainID() == GetChainID() && AMINO_ACID.GetSeqID() == GetLastAA()->GetSeqID() + 1,
        "The given amino acid: " + AMINO_ACID.GetIdentification() + " cannot be added to this SSE " +
        GetIdentification()
      );

      // make a copy of this amino acid
      util::ShPtr< biol::AABase> new_aa( AMINO_ACID.Clone());

      // initialize transformations
      math::TransformationMatrix3D transform;
      math::TransformationMatrix3D relative_residue_transform;

      // if the SSE should be idealized
      if( IDEALIZE)
      {
        // idealize this SSE
        SetToIdealConformationInPlace();

        // create a transformation to store the relative residue transformation
        relative_residue_transform
        (
          ( GetSize() * m_SSType->GetRiseInZPerResidue() / 2.0) * ( *coord::GetAxes().e_Z)
        );

        // apply the correct rotation
        relative_residue_transform( coord::GetAxes().e_Z, GetSize() * m_SSType->GetAnglePerTurn());

        // now to transform to the original coordinate
        transform = this->GetOrientation();
      }
      // do not idealize SSE and determine transformation from last residue
      else
      {
        // get the last aa
        util::ShPtr< biol::AABase> last_aa( GetData().LastElement());

        // make a hardcopy and set to ideal at origin
        util::ShPtr< biol::AABase> copy_aa( last_aa.HardCopy());
        copy_aa->SetToIdealConformation( m_SSType, math::TransformationMatrix3D());

        // get the transformation matrix to superimpose the aa's
        transform = quality::RMSD::SuperimposeCoordinates( last_aa->GetAtomCoordinates(), copy_aa->GetAtomCoordinates());

        // apply the incremental residue transformation
        relative_residue_transform( m_SSType->GetTransformationMatrixForResidues());
      }

      // now set to ideal coordinates for the calculated transformation
      new_aa->SetToIdealConformation( m_SSType, relative_residue_transform);

      // now to transform to the original coordinate
      new_aa->Transform( transform);

      // prepend the amino acid into the sequence
      biol::AASequence::PushBack( *new_aa);

      // set geometry
      SetGeometry();
    }

    //! @brief calculates ideal coordinates for the amino acids in the given sequence and prepends them to the sequence
    //! @param AA_SEQUENCE sequence to be prepended
    //! @param IDEALIZE bool whether to idealize the SSE prior to prepending
    void SSE::PrependSequence( const AASequence &AA_SEQUENCE, const bool IDEALIZE)
    {
      // if empty sequence is passed then do nothing
      if( AA_SEQUENCE.GetSize() == 0)
      {
        return;
      }

      // check the chain ids match, as well as seq id ares in correct order and the given sequence is continous
      BCL_Assert
      (
        AA_SEQUENCE.DoesPrecede( *this),
        "The given sequence: " + AA_SEQUENCE.GetSequenceIdentification() + " cannot be prepended to this SSE " +
        GetIdentification()
      );

      // initialize the transformation matrix
      math::TransformationMatrix3D transform;

      // if the SSE should be idealized
      if( IDEALIZE)
      {
        // idealize this SSE
        SetToIdealConformationInPlace();

        // initialize the transformation matrix to the body of the fragment
        transform = this->GetOrientation();

        // find the translation
        transform( -1.0 * GetLength() / 2.0 * GetAxis( coord::GetAxes().e_Z));
      }
      // do not idealize SSE and determine transformation from last residue
      else
      {
        // get the last aa
        util::ShPtr< biol::AABase> first_aa( GetData().FirstElement());

        // make a hardcopy and set to ideal at origin
        util::ShPtr< biol::AABase> copy_aa( first_aa.HardCopy());
        copy_aa->SetToIdealConformation( m_SSType, math::TransformationMatrix3D());

        // get the transformation matrix to superimpose the aa's
        transform = quality::RMSD::SuperimposeCoordinates( first_aa->GetAtomCoordinates(), copy_aa->GetAtomCoordinates());
      }

      // make a copy for the transformation for going from one residue to the next and invert it
      math::TransformationMatrix3D transform_for_prev_residue( m_SSType->GetTransformationMatrixForResidues());
      transform_for_prev_residue.Invert();

      // initialize a transformation for the residue to build up the cumulative transformation
      math::TransformationMatrix3D cumulative_transform_for_residue;

      // iterate over residues to add in reverse order
      for
      (
        biol::AASequence::const_reverse_iterator aa_itr( AA_SEQUENCE.ReverseBegin()),
          aa_itr_end( AA_SEQUENCE.ReverseEnd());
        aa_itr != aa_itr_end; ++aa_itr
      )
      {
        // make a copy of this amino acid
        util::ShPtr< biol::AABase> new_aa( ( *aa_itr)->Clone());

        // update the cumulative residue transformation
        cumulative_transform_for_residue( transform_for_prev_residue);

        // now set to ideal coordinates for the calculated transformation at the origin
        new_aa->SetToIdealConformation( m_SSType, cumulative_transform_for_residue);

        // now to transform to the original coordinate
        new_aa->Transform( transform);

        // prepend the amino acid into the sequence
        biol::AASequence::PushFront( *new_aa);
      }

      // set geometry
      SetGeometry();
    }

    //! @brief calculates ideal coordinates for the amino acids in the given sequence and appends them to the sequence
    //! @param AA_SEQUENCE sequence acid to be appended
    //! @param IDEALIZE bool whether to idealize the SSE prior to appending
    void SSE::AppendSequence( const AASequence &AA_SEQUENCE, const bool IDEALIZE)
    {
      // if empty sequence is passed then do nothing
      if( AA_SEQUENCE.GetSize() == 0)
      {
        return;
      }

      // if already peptide bonded, don't fix position of AA_SEQUENCE
      if( !IDEALIZE && biol::AABase::AreAminoAcidsPeptideBonded( *GetLastAA(), *AA_SEQUENCE.GetFirstAA(), true))
      {
        biol::AASequence::AppendSequence( AA_SEQUENCE);
      }
      else
      {
        // check the chain ids match, as well as seq id ares in correct order and the given sequence is continous
        BCL_Assert
        (
          DoesPrecede( AA_SEQUENCE),
          "The given sequence: " + AA_SEQUENCE.GetSequenceIdentification() + " cannot be appended to this SSE " +
          GetIdentification()
        );

        // create a transformation to store the relative residue transformation
        math::TransformationMatrix3D relative_residue_transform;
        math::TransformationMatrix3D transform;

        // if the SSE should be idealized
        if( IDEALIZE)
        {
          // idealize this SSE
          SetToIdealConformationInPlace();

          relative_residue_transform
          (
            ( ( GetSize() - 1) * m_SSType->GetRiseInZPerResidue() - GetLength() / 2.0) * ( *coord::GetAxes().e_Z)
          );

          // apply the correct rotation and overall transformation
          relative_residue_transform( coord::GetAxes().e_Z, ( GetSize() - 1) * m_SSType->GetAnglePerTurn());
          transform = this->GetOrientation();
        }
        // do not idealize SSE and determine transformation from last residue
        else
        {
          // get the last aa
          const biol::AABase &last_aa( *GetData().LastElement());

          // make a hardcopy and set to ideal at origin
          util::ShPtr< biol::AABase> copy_aa( last_aa.Clone());
          copy_aa->SetToIdealConformation( m_SSType, math::TransformationMatrix3D());

          // get the transformation matrix to superimpose the aa's
          transform = quality::RMSD::SuperimposeCoordinates
              (
                last_aa.GetAtomCoordinates( biol::GetAtomTypes().GetBackBoneAtomTypes()),
                copy_aa->GetAtomCoordinates( biol::GetAtomTypes().GetBackBoneAtomTypes())
              );
        }

        // iterate over residues to add
        for
        (
          biol::AASequence::const_iterator aa_itr( AA_SEQUENCE.Begin()), aa_itr_end( AA_SEQUENCE.End());
          aa_itr != aa_itr_end; ++aa_itr
        )
        {
          // make a copy of this amino acid
          util::ShPtr< biol::AABase> new_aa( ( *aa_itr)->Clone());

          // apply the incremental residue transformation
          relative_residue_transform( m_SSType->GetTransformationMatrixForResidues());

          // now set to ideal coordinates for the calculated transformation
          new_aa->SetToIdealConformation( m_SSType, relative_residue_transform);

          // now to transform to the original coordinate
          new_aa->Transform( transform);

          // prepend the amino acid into the sequence
          biol::AASequence::PushBack( *new_aa);
        }
      }

      // set geometry
      SetGeometry();
    }

    //! @brief fits SSE using this sequence and the passed SSE as a template
    //! @param SSE_TEMPLATE SSE to be used as a template
    //! @return fit SSE
    void SSE::FitToSSE( const SSE &SSE_TEMPLATE)
    {
      // if the SSE's are different types
      if( m_SSType != SSE_TEMPLATE.GetType())
      {
        // warn the user and switch this type
        BCL_MessageCrt( "Fitting SSE to another type!");
        m_SSType = SSE_TEMPLATE.GetType();
      }

      // fit the sequence
      biol::AASequenceFactory::FitSequence
      (
        *this,
        biol::AASequencePhiPsi( SSE_TEMPLATE),
        m_SSType
      );

      // update the geometry
      SetGeometry();
    }

    //! @brief get the longest contiguous predicted TM-segment
    bool SSE::IsPredictedTransmembrane() const
    {
      size_t length_current_seg( 0);
      auto method( sspred::GetMethods().GetCommandLineMethods());
      for( auto itr( this->GetData().Begin()), itr_end( this->GetData().End()); itr != itr_end; ++itr)
      {
        bool was_in_mem( false);
        for( auto itr_method( method.Begin()), itr_method_end( method.End()); itr_method != itr_method_end; ++itr_method)
        {
          auto method_pred( ( *itr)->GetSSPrediction( *itr_method));
          if( method_pred.IsDefined() && method_pred->GetOneStateTMPrediction() == biol::GetEnvironmentTypes().e_MembraneCore)
          {
            was_in_mem = true;
            break;
          }
        }
        if( !was_in_mem)
        {
          length_current_seg = 0;
        }
        else
        {
          ++length_current_seg;
          if
          (
            ( length_current_seg > size_t( 16) && m_SSType == biol::GetSSTypes().HELIX)
            || ( length_current_seg > size_t( 5) && m_SSType == biol::GetSSTypes().STRAND)
          )
          {
            return true;
          }
        }
      }
      return false;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator = SSE : makes HardCopy of the m_Data through the BaseClass operator =
    //! @param SSE_RHS SSE to be assigned to
    //! @return this SSE after assignment
    SSE &SSE::operator =( const SSE &SSE_RHS)
    {
      biol::AASequence::operator =( SSE_RHS);
      m_SSType = SSE_RHS.m_SSType;
      m_Orientation = SSE_RHS.m_Orientation;
      m_XExtent = SSE_RHS.m_XExtent;
      m_YExtent = SSE_RHS.m_YExtent;
      m_Fragments = SSE_RHS.m_Fragments;

      // emit coordinate change signal
      m_GeometryCoordinateChangeSignal.Emit( *this);
      m_CoordinateChangeSignal.Emit( *this);

      // retrun this
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SSE::Read( std::istream &ISTREAM)
    {
      // read base class
      biol::AASequence::Read( ISTREAM);

      // read the members
      io::Serialize::Read( m_SSType, ISTREAM);
      io::Serialize::Read( m_Orientation, ISTREAM);
      io::Serialize::Read( m_XExtent, ISTREAM);
      io::Serialize::Read( m_YExtent, ISTREAM);
      io::Serialize::Read( m_Fragments, ISTREAM);

      // emit coordinate change signal
      m_GeometryCoordinateChangeSignal.Emit( *this);
      m_CoordinateChangeSignal.Emit( *this);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SSE::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write base class
      biol::AASequence::Write( OSTREAM, INDENT) << '\n';

      // write the members
      io::Serialize::Write( m_SSType, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Orientation, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_XExtent, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_YExtent, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Fragments, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief sets the fragment geometries
    //! @param FRAGMENT_LENGTH fragment length to use
    void SSE::SetFragmentGeometries( size_t FRAGMENT_LENGTH)
    {
      // set the default fragment length (if no FRAGMENT_LENGTH was given)
      if( FRAGMENT_LENGTH == 0)
      {
        FRAGMENT_LENGTH = m_SSType->GetFragmentLength();
      }

      // if the sse is not the right type or the sse is the minimal size or smaller, return
      if( !m_SSType->IsStructured() || GetSize() <= FRAGMENT_LENGTH)
      {
        return;
      }

      // clear the fragment list
      m_Fragments.Reset();

      // calculate the number of geometries
      const size_t number_geometries( GetSize() - FRAGMENT_LENGTH + 1);

      // iterate over the number of geometries
      for( size_t geometry_ctr( 0); geometry_ctr < number_geometries; ++geometry_ctr)
      {
        // get the subsequence that corresponds to this fragment
        biol::AASequence sub_sequence
        (
          SubSequence
          (
            geometry_ctr, // starting position
            FRAGMENT_LENGTH // ending position
          )
        );

        // get central amino acid
        const int center( sub_sequence.GetFirstMember()->GetSeqID() + ( FRAGMENT_LENGTH - 1) / 2);

        // pushback the sse fragment into the fragments vector
        m_Fragments.PushBack
        (
          util::ShPtr< SSEGeometryInterface>
          (
            new SSEGeometry( SSE( sub_sequence, m_SSType), center)
          )
        );
      }

      // emit coordinate change signal
      m_GeometryCoordinateChangeSignal.Emit( *this);
      m_CoordinateChangeSignal.Emit( *this);
    }

    //! @brief sets the extents for the geometry
    void SSE::SetExtents()
    {
      // set the extents
      m_XExtent = m_SSType->GetRadialExtent();
      m_YExtent = m_SSType->GetRadialExtent();
    }

    //! @brief transforms only the geometries while leaving the AASequence coords in place
    //! @param TRANSFORMATION_MATRIX_3D transformation to apply to each geometry
    void SSE::TransformGeometries( const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D)
    {
      // transform orientation
      m_Orientation( TRANSFORMATION_MATRIX_3D);

      // transform fragment geometries
      for
      (
        util::ShPtrVector< SSEGeometryInterface>::iterator fragment_itr( m_Fragments.Begin()),
          fragment_itr_end( m_Fragments.End());
        fragment_itr != fragment_itr_end;
        ++fragment_itr
      )
      {
        ( *fragment_itr)->Transform( TRANSFORMATION_MATRIX_3D);
      }

      // emit coordinate change signal
      m_GeometryCoordinateChangeSignal.Emit( *this);
      m_CoordinateChangeSignal.Emit( *this);
    }

    //! @brief boolean operator SSE_LHS == SSE_RHS
    //! @param SSE_LHS first SSE
    //! @param SSE_RHS second SSE
    //! @return whether SSE_LHS is equal to SSE_RHS
    bool operator ==( const SSE &SSE_LHS, const SSE &SSE_RHS)
    {
      return
      (
        SSE_LHS.GetChainID() == SSE_RHS.GetChainID() &&
        SSE_LHS.GetType()    == SSE_RHS.GetType()    &&
        SSE_LHS.Sequence()   == SSE_RHS.Sequence()   &&
        SSE_LHS.GetFirstAA()->GetSeqID() == SSE_RHS.GetFirstAA()->GetSeqID() &&
        SSE_LHS.GetLastAA()->GetSeqID() == SSE_RHS.GetLastAA()->GetSeqID()
      );
    }

  } // namespace assemble
} // namespace bcl

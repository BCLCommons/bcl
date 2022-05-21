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
#include "restraint/bcl_restraint_cone_model.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_collector_aa_type.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "quality/bcl_quality_rmsd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ConeModel::s_Instance
    (
      GetObjectInstances().AddInstance( new ConeModel())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ConeModel::ConeModel()
    {
    }

    //! @brief Clone function
    //! @return pointer to new ConeModel
    ConeModel *ConeModel::Clone() const
    {
      return new ConeModel( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ConeModel::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

      //! @brief calculates the maximum angle (SL->CB->SL) observed for an ensemble of models
      //! @param ENSEMBLE ensemble for which the maximum SL->CB->SL angle will be determined
      //! @return the maximum SL->CB->SL angle observed in the provided ensemble
      double ConeModel::SLCBSLMaxAngle( const assemble::ProteinEnsemble &ENSEMBLE)
      {

        // will hold the maximum SL->SB->SL angle seen for this ensemble
        double sl_cb_sl_max_angle( 0);

        // iterate through the ensemble
        for
        (
          assemble::ProteinEnsemble::const_iterator model_itr_a( ENSEMBLE.Begin()), model_itr_end( ENSEMBLE.End());
          model_itr_a != model_itr_end; ++model_itr_a
        )
        {
          // get the spin label residue from model_a
          const util::SiPtr< const biol::AABase> sl_a( GetSpinLabelResidue( **model_itr_a));

          // get the atom coordinates of the spin label side chain
          const util::SiPtrVector< const linal::Vector3D> sl_a_coords
          (
            sl_a->GetAtomCoordinates( GetSuperimposeAtomTypes())
          );

          // get the unpaired electron coordinates from spin label a
          const linal::Vector3D unpaired_electron_a_coords( GetUnpairedElectronCoordinates( *sl_a));

          // iterate through the ensemble a second time
          assemble::ProteinEnsemble::const_iterator model_itr_b( model_itr_a);
          ++model_itr_b;
          for( ; model_itr_b != model_itr_end; ++model_itr_b)
          {
            // get the spin label coordinates from spin label b copy that has been superimposed on spin label a
            const linal::Vector3D unpaired_electron_b_coords( GetSuperimposedUnpairedElectronCoordinates( **model_itr_b, sl_a_coords));

            // calculate the angle between the CB of sl_a, and the electron in sl_a and sl_b_copy
            const double proj_angle
            (
              linal::ProjAngle
              (
                sl_a->GetAtom( biol::GetAtomTypes().CB).GetCoordinates(),
                unpaired_electron_a_coords,
                unpaired_electron_b_coords
              )
            );

            BCL_MessageDbg( "current SL->CB-SL angle is " + util::Format()( proj_angle));

            // true if the current SL->CB->SL angle is larger than max angle seen so far
            if( proj_angle > sl_cb_sl_max_angle)
            {
              // set the max angle to the current projection angle
              sl_cb_sl_max_angle = proj_angle;
            }
          }
        }

        // return the maximum SL->CB->SL angle observed in the ensemble
        return sl_cb_sl_max_angle;
      }

      //! @brief calculates the angle between the effective spin label position, the CB, and the CA
      //!        SLeffective->CB->CA
      //! @param ENSEMBLE the ensemble of protein models from which the SLeffective->CB->CA
      //! @return double which is the calculated SLeffective->CB->CA angle
      double ConeModel::SLeffectiveCBCAAngle( const assemble::ProteinEnsemble &ENSEMBLE)
      {
        storage::Pair< linal::Vector3D, util::SiPtr< const biol::AABase> >
          slef_tmpltaa( GetSLEffectivePositionAndTemplateResidue( ENSEMBLE));

        // now calculate the SLeffective->CB->CA angle using the sl_a coordinates for CB and CA
        const double proj_angle
        (
          linal::ProjAngle
          (
            slef_tmpltaa.Second()->GetAtom( biol::GetAtomTypes().CB).GetCoordinates(),
            slef_tmpltaa.Second()->GetAtom( biol::GetAtomTypes().CA).GetCoordinates(),
            slef_tmpltaa.First()
          )
        );

        BCL_MessageDbg( "SLeffective->CB->CA angle is " + util::Format()( proj_angle));

        // return the calculated SLeffective->CB->CA angle
        return proj_angle;
      }

      //! @brief calculates the distance between the effective spin label position and the CB for an ensemble
      //! @param ENSEMBLE the protein ensemble for which the SL effective position will be calculated
      //! @return double which is the distance between the CB and SLeffective position
      double ConeModel::SLeffectiveCBDistance( const assemble::ProteinEnsemble &ENSEMBLE)
      {
        const storage::Pair< linal::Vector3D, util::SiPtr< const biol::AABase> > sleff_tmpltaa
        (
          GetSLEffectivePositionAndTemplateResidue( ENSEMBLE)
        );

        const double distance
        (
          linal::Distance
          (
            sleff_tmpltaa.First(), sleff_tmpltaa.Second()->GetAtom( biol::GetAtomTypes().CB).GetCoordinates()
          )
        );

        BCL_MessageDbg( "current SLeffective->CB distance is " + util::Format()( distance));

        return distance;
      }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ConeModel::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ConeModel::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

      //! @brief calculates the position of the unpaired electron in a spin label residue
      //!        the midpoint of the bond between the oxygen and nitrogen in the nitroxide moiety in the ring
      //! @param RESIDUE the residue where the position of the spin label will be calculated
      //! @return vector 3d which has the coordinates of the spin label (i.e. unpaired electron)
      linal::Vector3D ConeModel::GetUnpairedElectronCoordinates( const biol::AABase &RESIDUE)
      {
        // get the coordinates for the oxygen and nitrogen in the nitroxide moiety
        const util::SiPtrVector< const linal::Vector3D> o_coords
        (
          RESIDUE.GetAtomCoordinates( storage::Set< biol::AtomType>( biol::GetAtomTypes().O1))
        );
        BCL_Assert
        (
          o_coords.GetSize() == 1, util::Format()( o_coords.GetSize()) + " O1 coords found in " +
          RESIDUE.GetIdentification() + "\n" + util::Format()( RESIDUE) + "\nBe sure you are using AAComplete"
        );
        const util::SiPtrVector< const linal::Vector3D> n_coords
        (
          RESIDUE.GetAtomCoordinates( storage::Set< biol::AtomType>( biol::GetAtomTypes().N1))
        );

        BCL_Assert
        (
          n_coords.GetSize() == 1, util::Format()( n_coords.GetSize()) + " N1 coords found in " +
          RESIDUE.GetIdentification() + "\n" + util::Format()( RESIDUE) + "\nBe sure you are using AAComplete"
        );

        // calculate the position of the unpaired electron
        const linal::Vector3D ave_pos( ( *o_coords.FirstElement() + *n_coords.FirstElement()) / 2.0);

        // return SL coordinates
        return ave_pos;
      }

      //! @brief gives the single spin label residue that is within a protein model - asserts if there is more than one
      //! @param MODEL the model which has one spin label residue which will be given
      //! @return SiPtr to the spin label residue that is within the protein model
      util::SiPtr< const biol::AABase> ConeModel::GetSpinLabelResidue
      (
        const assemble::ProteinModel &MODEL
      )
      {
        // will collect the spin label residues from the protein models
        const assemble::CollectorAAType sl_collector( storage::Set< biol::AAType>( biol::GetAATypes().R1A));

        // get the spin label residue from model
        const util::SiPtrList< const biol::AABase> sl_list_a( sl_collector.Collect( MODEL.GetAminoAcids()));

        // make sure exactly one spin label side residue was found
        BCL_Assert( sl_list_a.GetSize() == 1, util::Format()( sl_list_a.GetSize()) + " spin labels found");

        // get the one spin label residue out of the list
        const util::SiPtr< const biol::AABase> sl_a( sl_list_a.FirstElement());

        // return the coordinates of the unpaired electron
        return sl_a;
      }

      //! @brief gives the atom types that are used to superimpose spin label residues when the effective spin label
      //!        position is being determined
      //! @return set which has the atom types used for superimposition
      const storage::Set< biol::AtomType> &ConeModel::GetSuperimposeAtomTypes()
      {
        static const storage::Set< biol::AtomType> s_types
        (
          storage::Set< biol::AtomType>::Create
          (
            biol::GetAtomTypes().CA, biol::GetAtomTypes().CB, biol::GetAtomTypes().N, biol::GetAtomTypes().C,
            biol::GetAtomTypes().O
          )
        );

        return s_types;
      }

      //! @brief gets a spin label residue from a protein model and superimposes it onto the provided coordinates
      //!        and returns the resulting position of the unpaired electron in that spin label
      //! @param MODEL the model from which the spin label residue is going to be gotten
      //! @param COORDS_TO_SUPERIMPOSE_ON the coordinates to which the spin label residue will be superimposed
      //! @return Vector3D which has the coordinates of the unpaired electron in the superimposed SL residue
      linal::Vector3D ConeModel::GetSuperimposedUnpairedElectronCoordinates
      (
        const assemble::ProteinModel &MODEL, const util::SiPtrVector< const linal::Vector3D> &COORDS_TO_SUPERIMPOSE_ON
      )
      {
        // get the spin label residue from model
        const util::SiPtr< const biol::AABase> sl( GetSpinLabelResidue( MODEL));

        // get the atom coordinates of the spin label b side chain
        const util::SiPtrVector< const linal::Vector3D> sl_coords( sl->GetAtomCoordinates( GetSuperimposeAtomTypes()));

        // get the transformation matrix necessary to superimpose the CA, CB, N, and C atoms of the two sl residues
        // transformation matrix that will superimpose sl_b_coords onto sl_a_coords
        const math::TransformationMatrix3D transform( quality::RMSD::SuperimposeCoordinates( COORDS_TO_SUPERIMPOSE_ON, sl_coords));

        // copy spin label b and superimpose sl_b onto sl_a
        util::ShPtr< biol::AABase> sl_copy( sl->Clone());
        sl_copy->Transform( transform);

        // get the spin label coordinates from spin label b copy that has been superimposed on spin label a
        const linal::Vector3D unpaired_electron_coords( GetUnpairedElectronCoordinates( *sl_copy));

        return unpaired_electron_coords;
      }

      //! @brief gives the coordinate of the effective spin label (unpaired electron) position and the amino acid
      //!        used as the basis for superimposing all the spin labels into the same space
      //! @param ENSEMBLE the ensemble for which the effective spin label position will be calculated
      //! @return Pair which has the coordinate of the effective spin label and the spin label used for superimposition
      storage::Pair< linal::Vector3D, util::SiPtr< const biol::AABase> >
      ConeModel::GetSLEffectivePositionAndTemplateResidue( const assemble::ProteinEnsemble &ENSEMBLE)
      {
        linal::Vector3D sl_effective_pos( 0, 0, 0);

        util::SiPtr< const biol::AABase> sl_a( GetSpinLabelResidue( **ENSEMBLE.Begin()));

        // get the atom coordinates of the spin label side chain
        const util::SiPtrVector< const linal::Vector3D> sl_a_coords
        (
          sl_a->GetAtomCoordinates( GetSuperimposeAtomTypes())
        );

        sl_effective_pos += GetUnpairedElectronCoordinates( *sl_a);

        // iterate through the ensemble - starting with second protein, since first was already counted
        for
        (
          assemble::ProteinEnsemble::const_iterator model_itr( ++ENSEMBLE.Begin()), model_itr_end( ENSEMBLE.End());
          model_itr != model_itr_end; ++model_itr
        )
        {
          sl_effective_pos += GetSuperimposedUnpairedElectronCoordinates( **model_itr, sl_a_coords);
        }

        // get the average position of all the unpaired electrons, this is the effective SL position
        sl_effective_pos /= double( ENSEMBLE.GetSize());

        BCL_MessageDbg( "effective sl position is " + util::Format()( sl_effective_pos));

        return storage::Pair< linal::Vector3D, util::SiPtr< const biol::AABase> >( sl_effective_pos, sl_a);
      }

  } // namespace restraint

} // namespace bcl

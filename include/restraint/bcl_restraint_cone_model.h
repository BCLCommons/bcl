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

#ifndef BCL_RESTRAINT_CONE_MODEL_H_
#define BCL_RESTRAINT_CONE_MODEL_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "biol/bcl_biol.fwd.hh"
#include "linal/bcl_linal.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConeModel
    //! @brief has functionality for calculating parameters relevant to the spin label cone model
    //! @details Can be used to calculate cone model parameters.
    //!          See Alexander et al. Structure 2008 "De novo high-resolution protein structure determination from
    //!          sparse spin-labeling EPR data"
    //!          See Hirst et al. J. Struct. Bio. 2011 "RosettaEPR: An integrated tool for protein structure
    //!          determination from sparse EPR data"
    //!
    //! @see @link example_restraint_cone_model.cpp @endlink
    //! @author alexanns
    //! @date Mar 4, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ConeModel :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ConeModel();

      //! @brief Clone function
      //! @return pointer to new ConeModel
      ConeModel *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief calculates the maximum angle (SL->CB->SL) observed for an ensemble of models
      //! @param ENSEMBLE ensemble for which the maximum SL->CB->SL angle will be determined
      //! @return the maximum SL->CB->SL angle observed in the provided ensemble
      static double SLCBSLMaxAngle( const assemble::ProteinEnsemble &ENSEMBLE);

      //! @brief calculates the angle between the effective spin label position, the CB, and the CA
      //!        SLeffective->CB->CA
      //! @param ENSEMBLE the ensemble of protein models from which the SLeffective->CB->CA
      //! @return double which is the calculated SLeffective->CB->CA angle
      static double SLeffectiveCBCAAngle( const assemble::ProteinEnsemble &ENSEMBLE);

      //! @brief calculates the distance between the effective spin label position and the CB for an ensemble
      //! @param ENSEMBLE the protein ensemble for which the SL effective position will be calculated
      //! @return double which is the distance between the CB and SLeffective position
      static double SLeffectiveCBDistance( const assemble::ProteinEnsemble &ENSEMBLE);

    ///////////////
    // operators //
    ///////////////

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

    public:

      //! @brief calculates the position of the unpaired electron in a spin label residue
      //!        the midpoint of the bond between the oxygen and nitrogen in the nitroxide moiety in the ring
      //! @param RESIDUE the residue where the position of the spin label will be calculated
      //! @return vector 3d which has the coordinates of the spin label (i.e. unpaired electron)
      static linal::Vector3D GetUnpairedElectronCoordinates( const biol::AABase &RESIDUE);

      //! @brief gives the single spin label residue that is within a protein model - asserts if there is more than one
      //! @param MODEL the model which has one spin label residue which will be given
      //! @return SiPtr to the spin label residue that is within the protein model
      static util::SiPtr< const biol::AABase> GetSpinLabelResidue
      (
        const assemble::ProteinModel &MODEL
      );

    private:

      //! @brief gives the atom types that are used to superimpose spin label residues when the effective spin label
      //!        position is being determined
      //! @return set which has the atom types used for superimposition
      static const storage::Set< biol::AtomType> &GetSuperimposeAtomTypes();

      //! @brief gets a spin label residue from a protein model and superimposes it onto the provided coordinates
      //!        and returns the resulting position of the unpaired electron in that spin label
      //! @param MODEL the model from which the spin label residue is going to be gotten
      //! @param COORDS_TO_SUPERIMPOSE_ON the coordinates to which the spin label residue will be superimposed
      //! @return Vector3D which has the coordinates of the unpaired electron in the superimposed SL residue
      static linal::Vector3D GetSuperimposedUnpairedElectronCoordinates
      (
        const assemble::ProteinModel &MODEL, const util::SiPtrVector< const linal::Vector3D> &COORDS_TO_SUPERIMPOSE_ON
      );

      //! @brief gives the coordinate of the effective spin label (unpaired electron) position and the amino acid
      //!        used as the basis for superimposing all the spin labels into the same space
      //! @param ENSEMBLE the ensemble for which the effective spin label position will be calculated
      //! @return Pair which has the coordinate of the effective spin label and the spin label used for superimposition
      static storage::Pair< linal::Vector3D, util::SiPtr< const biol::AABase> >
      GetSLEffectivePositionAndTemplateResidue( const assemble::ProteinEnsemble &ENSEMBLE);

    }; // class ConeModel

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_CONE_MODEL_H_ 

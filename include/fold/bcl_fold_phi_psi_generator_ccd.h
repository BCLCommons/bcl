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

#ifndef BCL_FOLD_PHI_PSI_GENERATOR_CCD_H_
#define BCL_FOLD_PHI_PSI_GENERATOR_CCD_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_fold_mutation_residue.h"
#include "biol/bcl_biol_aa_sequence_flexibility.h"
#include "coord/bcl_coord_cyclic_coordinate_descent.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "random/bcl_random_distribution_interface.h"
#include "storage/bcl_storage_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PhiPsiGeneratorCCD
    //! @brief is for generating phi and psi angles such that the phi or psi angle it generates will minimize the sum
    //! square distance deviation between a list of target a moving points.
    //! @details  It mutates either phi or psi and gives an undefined double for the other. Mutating phi or psi is
    //! randomly determined. The determined optimal phi or psi value itself is not used in the end, but a randomly
    //! chosen fraction of it is.
    //!
    //! @see @link example_fold_phi_psi_generator_ccd.cpp @endlink
    //! @author alexanns
    //! @date Sep 6, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PhiPsiGeneratorCCD :
      public math::FunctionInterfaceSerializable< MutationResidue, storage::VectorND< 2, double> >
    {

    private:

    //////////
    // data //
    //////////

      //! the list of target and moving points whose distances will be minimized
      storage::List< coord::CyclicCoordinateDescent::TargetAndMovingPointPair> m_TargetAndMovingPoints;

      //! the random number generator that should be used
      const random::DistributionInterface &m_RandomNumberGenerator;

      //! indicates if the direction the rotation needs to occur in
      biol::AASequenceFlexibility::DirectionEnum m_Direction;

      //! range for a random fraction to multiply with the suggested rotation angle
      math::Range< double> m_RandomFraction;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      PhiPsiGeneratorCCD();

      //! @brief constructor taking member variable parameters
      //! @param TARGET_AND_MOVING_POINTS the list of target and moving points whose distances will be minimized
      //! @param RANDOM_NUMBER_GENERATOR the random number generator that should be used
      //! @param DIRECTION indicates if the rotation needs to occur in the C to N direction
      //! @param RANDOM_FRACTION_RANGE range a random fraction is drawn from and multiplied with suggested rotation
      PhiPsiGeneratorCCD
      (
        const storage::List< coord::CyclicCoordinateDescent::TargetAndMovingPointPair> &TARGET_AND_MOVING_POINTS,
        const random::DistributionInterface &RANDOM_NUMBER_GENERATOR,
        const biol::AASequenceFlexibility::SequenceDirection &DIRECTION,
        const math::Range< double> &RANDOM_FRACTION_RANGE
      );

      //! @brief Clone function
      //! @return pointer to new PhiPsiGeneratorCCD
      PhiPsiGeneratorCCD *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator taking MutationResidue and returning storage::VectorND< 2, double> with phi psi, respectively
      //! @param MUTATION_RESIDUE MutationResidue whose phi or psi will be changed
      //! @return phi or psi value, respectively, which will minimize the distance between target and moving points.
      //!         the other value will be undefined
      storage::VectorND< 2, double> operator()( const MutationResidue &MUTATION_RESIDUE) const;

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

    }; // class PhiPsiGeneratorCCD

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_PHI_PSI_GENERATOR_CCD_H_ 

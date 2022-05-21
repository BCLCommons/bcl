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

#ifndef BCL_FOLD_MUTATE_AA_SET_PHI_H_
#define BCL_FOLD_MUTATE_AA_SET_PHI_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_atom.h"
#include "math/bcl_math_mutate_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateAASetPhi
    //! @brief is for setting the phi angle of an amino acid to a desired angle.
    //! details  It rotates around the N->CA bond such that the N and CA atoms are stationary but all other atoms
    //! c-terminal are moved as necessary so that the phi angle is set as needed.
    //!
    //! @see @link example_fold_mutate_aa_set_phi.cpp @endlink
    //! @author alexanns
    //! @date Sep 4, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateAASetPhi :
      public math::MutateInterface< biol::AABase>
    {

    private:

    //////////
    // data //
    //////////

      //! the C atom of the residue previous to the residue whose phi will be changed
      biol::Atom m_PreviousC;

      //! the desired phi angle value in radians
      double m_DesiredPhi;

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
      MutateAASetPhi();

      //! @brief constructor taking member variables
      //! @param PREVIOUS_C the C atom of the residue previous to the residue whose phi will be changed
      //! @param DESIRED_PHI the desired phi angle value in radians
      MutateAASetPhi( const biol::Atom &PREVIOUS_C, const double DESIRED_PHI);

      //! @brief Clone function
      //! @return pointer to new MutateAASetPhi
      MutateAASetPhi *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief GetTransformationMatrix gives the transformation matrix that is necessary to set phi of a residue
      //! @param RESIDUE the residue for which the transformation matrix will be calculated
      //! @return the transformation matrix that is necessary to set phi of "RESIDUE" to "m_DesiredPhi"
      math::TransformationMatrix3D GetTransformationMatrix( const biol::AABase &RESIDUE) const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator taking an ARGUMENT and returning a mutated object of t_ArgumentType
      //! @param RESIDUE Argument of interest
      //! @return MutateResult that results from mutating to the argument
      math::MutateResult< biol::AABase> operator()( const biol::AABase &RESIDUE) const;

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

    protected:

      //! @brief GetMutateAAPhi provides the MutateAAPhi object necesssary to set phi as desired
      //! @param RESIDUE is the residue whose phi is going to be set
      //! @return the MutateAAPhi object necesssary to set phi as desired
      MutateAAPhi GetMutateAAPhi( const biol::AABase &RESIDUE) const;

    }; // class MutateAASetPhi

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_AA_SET_PHI_H_ 

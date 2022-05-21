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

#ifndef BCL_FOLD_MUTATE_AA_SET_PSI_H_
#define BCL_FOLD_MUTATE_AA_SET_PSI_H_

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
    //! @class MutateAASetPsi
    //! @brief is for setting the psi angle of an amino acid to a desired angle. The rotation is
    //! performed around the CA->C bond vector, so the CA and C positions are fixed. The oxygen atom attached
    //! to the C atom is rotated accordingly.
    //!
    //! @see @link example_fold_mutate_aa_set_psi.cpp @endlink
    //! @author alexanns
    //! @date Sep 4, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateAASetPsi :
      public math::MutateInterface< biol::AABase>
    {

    private:

    //////////
    // data //
    //////////

      //! the N atom of the residue following the residue whose psi will be changed
      biol::Atom m_FollowingN;

      //! the desired psi angle value in radians
      double m_DesiredPsi;

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
      MutateAASetPsi();

      //! @brief constructor taking member variables
      //! @param FOLLOWING_N the N atom of the residue following the residue whose psi will be changed
      //! @param DESIRED_PSI the desired psi angle value in radians
      MutateAASetPsi( const biol::Atom &FOLLOWING_N, const double DESIRED_PSI);

      //! @brief Clone function
      //! @return pointer to new MutateAASetPsi
      MutateAASetPsi *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief GetTransformationMatrix gives the transformation matrix that is necessary to set psi of a residue
      //! @param RESIDUE the residue for which the transformation matrix will be calculated
      //! @return the transformation matrix that is necessary to set psi of "RESIDUE" to "m_DesiredPsi"
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

      //! @brief GetMutateAAPsi provides the MutateAAPsi object necesssary to set psi as desired
      //! @param RESIDUE is the residue whose psi is going to be set
      //! @return the MutateAAPsi object necesssary to set psi as desired
      MutateAAPsi GetMutateAAPsi( const biol::AABase &RESIDUE) const;

    }; // class MutateAASetPsi

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_AA_SET_PSI_H_ 

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

#ifndef BCL_CHEMISTRY_BOND_DIHEDRAL_STERIC_WEIGHT_H_
#define BCL_CHEMISTRY_BOND_DIHEDRAL_STERIC_WEIGHT_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_complete.h"
#include "bcl_chemistry_atom_vector.h"
#include "linal/bcl_linal_vector.h"
#include "sdf/bcl_sdf_molecule_reading_pref.h"
#include "storage/bcl_storage_triplet.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BondDihedralStericWeight
    //! @brief Standardizes FragmentComplete object
    //!
    //! @see @link example_chemistry_bond_dihedral_steric_weight.cpp @endlink
    //! @author mendenjl
    //! @date Oct 04, 2017
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API BondDihedralStericWeight :
      public util::ObjectInterface
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new BondDihedralStericWeight
      BondDihedralStericWeight *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief Determines steric weight at each side of the bond separately
      //! @param CONFORMATION the conformations of interest
      //! @param N_SPHERES number of spheres to go out
      static storage::Vector< storage::Triplet< size_t, size_t, linal::Vector< double> > >
      CalculateDualSidedWeight
      (
        const FragmentComplete &CONFORMATION,
        const size_t N_SPHERES
      );

      //! @brief Determines steric weight at each side of the bond separately
      //! @param CONFORMATION the conformations of interest
      //! @param N_SPHERES number of spheres to go out
      static storage::Vector< graph::UndirectedEdge< double> >
      CalculateBondWeights
      (
        const FragmentComplete &CONFORMATION,
        const size_t N_SPHERES,
        const double DEPLETION_FACTOR
      );

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class BondDihedralStericWeight

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_BOND_DIHEDRAL_STERIC_WEIGHT_H_


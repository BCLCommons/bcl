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

#ifndef BCL_CHEMISTRY_ATOM_VDW_SCORE_H_
#define BCL_CHEMISTRY_ATOM_VDW_SCORE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_interface.h"
#include "math/bcl_math_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AtomVdwScore
    //! @brief Scores clashes between atoms of a given conformation
    //! @details Potential that evalutes clashes in a given molecule conformation
    //!
    //! @see @link example_chemistry_atom_vdw_score.cpp @endlink
    //! @author mendenjl
    //! @date Oct 25, 2019
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AtomVdwScore :
      public math::FunctionInterfaceSerializable< FragmentComplete, double>
    {

    private:

      //! Molecule dependent variables. These variables are cached because they are expensive to compute and the same
      //! molecule may have the clash score called on it hundreds of times in succession
      mutable storage::Vector< AtomType>      m_LastAtomTypes;
      mutable storage::Vector< sdf::BondInfo> m_LastBondInfo;

      mutable float                           m_MaxVdwRadius;
      mutable linal::Matrix< float>           m_VdwSumMaxDistances;
      mutable linal::Matrix< float>           m_PVdw;
      mutable linal::Matrix< float>           m_LJ;
      bool m_ClashOnly;
      mutable bool m_Display;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;
      static const util::SiPtr< const util::ObjectInterface> s_ClashInstance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      AtomVdwScore( const bool &CLASH_ONLY = false) :
        m_ClashOnly( CLASH_ONLY),
        m_Display( false)
      {
      }

      //! @brief Clone function
      //! @return pointer to new AtomVdwScore
      AtomVdwScore *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetAlias() const;

      void SetDisplay( const bool &DISP) const
      {
        m_Display = DISP;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief evaluate clashes for given atom pair
      //! @param MOLECULE molecule that needs to scored
      //! @return clash score for the given atom pair
      double operator()
      (
        const FragmentComplete &MOLECULE
      ) const;

      //! @brief evaluate clashes for given atoms
      //! @param MOLECULE molecule that needs to scored
      //! @return clash score for the given atoms
      double operator()
      (
        const ConformationInterface &MOLECULE
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    private:

      //! @brief Test whether this molecule is the same (constitutionally) as the molecule for which the state of this
      //!        class currently can handle
      bool IsMoleculeInfoCached( const ConformationInterface &CONF) const;

      //! @brief Update molecule change the molecule that this class will compute the clash score for
      //! @param MOL molecule of interest
      void UpdateMolecule( const ConformationInterface &CONF) const;

    }; // class AtomVdwScore

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_ATOM_VDW_SCORE_H_

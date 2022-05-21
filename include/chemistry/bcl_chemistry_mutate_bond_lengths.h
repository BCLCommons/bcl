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

#ifndef BCL_CHEMISTRY_MUTATE_BOND_LENGTHS_H_
#include "bcl_chemistry_bond_lengths.h"

#define BCL_CHEMISTRY_MUTATE_BOND_LENGTHS_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_interface.h"
#include "bcl_chemistry_mutate_bond_angles.h"
#include "bcl_chemistry_mutate_dihedral_bond.h"
#include "bcl_chemistry_mutate_dihedrals_interface.h"
#include "bcl_chemistry_rotamer_dihedral_bond_data.h"
#include "coord/bcl_coord_line_segment_3d.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_mutate_result.h"
#include "random/bcl_random_histogram_1d_distribution.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_own_ptr.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateBondLengths
    //! @brief This class adjusts bond lengths to equilibrium values by iteratively repositioning mobile atoms
    //! @details
    //!
    //! @see @link example_chemistry_mutate_bond_lengths.cpp @endlink
    //! @author brownbp1
    //! @date Oct 14, 2021
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateBondLengths :
      public math::MutateInterface< FragmentComplete>
    {

    private:

    ////////////////////
    // helper classes //
    ////////////////////

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class MolBondInfo
      //! @brief a container object to manage bond length accounting
      //!
      //! @author brownbp1
      //! @date Oct 14, 2021
      //!
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      struct MolBondInfo
      {
        //! atoms
        AtomVector< AtomComplete> m_Atoms;

        //! bondinfo
        storage::Vector< sdf::BondInfo> m_BondInfo;

        //! number of bonds
        size_t m_NumberBonds;

        //! equilibrium bond lengths
        storage::Vector< double> m_BondLengthsEq;

        //! current bond lengths
        storage::Vector< double> m_BondLengthsCurrent;

        //! bond tolerance per bond
        storage::Vector< double> m_BondLengthsTolerance;

        //! difference in equilibrium and current bond lengths accounting for tolerance
        storage::Vector< double> m_BondLengthsDiff;

        //! worst bond index by deviation from equilibrium
        size_t m_WorstBondIndex;

        //! @brief default construct
        MolBondInfo() :
          m_Atoms( AtomVector< AtomComplete>()),
          m_BondInfo( storage::Vector< sdf::BondInfo>()),
          m_NumberBonds( size_t()),
          m_BondLengthsEq( storage::Vector< double>()),
          m_BondLengthsCurrent( storage::Vector< double>()),
          m_BondLengthsTolerance( storage::Vector< double>()),
          m_BondLengthsDiff( storage::Vector< double>()),
          m_WorstBondIndex( size_t())
        {}

        //! @brief constructor with molecule and tolerance
        MolBondInfo
        (
          const FragmentComplete &MOL,
          const storage::Vector< double> &TOL = storage::Vector< double>()
        ) :
          m_Atoms( MOL.GetAtomVector()),
          m_BondInfo( MOL.GetBondInfo()),
          m_NumberBonds( MOL.GetBondInfo().GetSize()),
          m_BondLengthsEq( MOL.GetBondInfo().GetSize(), util::GetUndefinedDouble()),
          m_BondLengthsCurrent( MOL.GetBondInfo().GetSize(), util::GetUndefinedDouble()),
          m_BondLengthsTolerance( TOL),
          m_BondLengthsDiff( MOL.GetBondInfo().GetSize(), util::GetUndefinedDouble()),
          m_WorstBondIndex( util::GetUndefinedSize_t())
        {
          // set per-bond tolerances if not set properly already
          if( m_BondLengthsTolerance.IsEmpty() || m_BondLengthsTolerance.GetSize() != m_NumberBonds)
          {
            // default
            m_BondLengthsTolerance = storage::Vector< double>( m_NumberBonds, 0.1);
          }

          util::OwnPtr< int> blah( util::SiPtr< int>());

          // set values for each bond
          for( size_t i( 0); i < m_NumberBonds; ++i)
          {
            // get atom reference
            auto &atom_a( m_Atoms( m_BondInfo( i).GetAtomIndexLow()));
            auto &atom_b( m_Atoms( m_BondInfo( i).GetAtomIndexHigh()));

            // equilibrium bond length
            m_BondLengthsEq( i) = BondLengths::GetBondLength
            (
              atom_a.GetAtomType(),
              m_BondInfo( i).GetConfigurationalBondType()->GetBondData( ConfigurationalBondTypeData::e_BondOrderOrAromatic),
              atom_b.GetAtomType()
            );

            // current bond length
            linal::Vector3D bond_vec( atom_a.GetPosition() - atom_b.GetPosition());
            m_BondLengthsCurrent( i) = bond_vec.Norm();

            // difference from equilibrium accounting for tolerance
            m_BondLengthsDiff( i) =
                math::Absolute( m_BondLengthsCurrent( i) - m_BondLengthsEq( i)) > m_BondLengthsTolerance( i) ?
                    m_BondLengthsCurrent( i) - m_BondLengthsEq( i) :
                    0.0;
          }

          // get the worst bond index
          double greatest_deviation( 0.0);
          for( size_t i( 0); i < m_NumberBonds; ++i)
          {
            if( math::Absolute( m_BondLengthsDiff( i)) > greatest_deviation)
            {
              greatest_deviation = math::Absolute( m_BondLengthsDiff( i));
              m_WorstBondIndex = i;
            }
          }

          // if all of the bonds are good then set worst undefined
          if( !greatest_deviation)
          {
            m_WorstBondIndex = util::GetUndefinedSize_t();
          }
        } // end MolBondInfo constructor
      }; // end MolBondInfo

    //////////
    // data //
    //////////

      //! maximum number of iterations in which to adjust bond lengths
      size_t m_MaxCycles;

      //! bonds lengths +/- individual tolerance must be below this value for convergence
      double m_ConvergenceLimit;

      //! increment to use when adjusting bond lengths (Angstroms)
      double m_StepSize; // TODO consider making dynamic range

      //! atoms that are able to be moved during bond length adjustment
      mutable storage::Vector< size_t> m_MobileAtoms;

      //! whether the info object has been initialized
      mutable bool m_Initialized;

      //! info object for current molecule
      mutable MolBondInfo m_Info;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor taking the member variable parameters
      MutateBondLengths
      (
        const size_t &MAX_CYCLES = 5,
        const double &CONVERGENCE_LIMIT = 0.01,
        const double &STEP_SIZE = 0.01,
        const storage::Vector< size_t> &MOBILE_ATOMS = storage::Vector< size_t>()
      ) :
        m_MaxCycles( MAX_CYCLES),
        m_ConvergenceLimit( CONVERGENCE_LIMIT),
        m_StepSize( STEP_SIZE),
        m_MobileAtoms( MOBILE_ATOMS),
        m_Initialized( false),
        m_Info( MolBondInfo())
      {
      }

      //! @brief Clone function
      //! @return pointer to new MutateBondLengths
      MutateBondLengths *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetScheme() const;

      //! @brief get access to the stored molecule info
      //! @return the MolBondInfo object
      const MolBondInfo &GetInfo() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator taking a conformation and returns a mutated conformation
      //! @param MOLECULE conformation of interest
      //! @return MutateResult that results from mutating to the argument
      math::MutateResult< FragmentComplete> operator()
      (
        const FragmentComplete &MOLECULE
      ) const;

    ////////////////
    // operations //
    ////////////////

      //! @brief initialize info object with molecule
      void Initialize( const FragmentComplete &MOLECULE) const;

      // TODO: replace this with something better
      //! @brief setup this mutate to handle a new molecule
      storage::Vector< size_t> GetMobileAtoms
      (
        const FragmentComplete &CONF
      ) const;

      //! @brief perturb the molecule to update bond lengths
      //! @param INFO bond length accounting info of current molecule
      //! @param MOBILE_ATOMS atoms that can be moved during perturbation
      void UpdateBondLengths
      (
        MolBondInfo &INFO,
        const storage::Vector< size_t> MOBILE_ATOMS
      ) const;

      //! @brief update atom types of select atoms
      //! @param INFO bond types accounting info of current molecule
      //! @param MOBILE_ATOMS atoms whose bond type can be altered
      void UpdateAtomTypes
      (
        MolBondInfo &INFO,
        const storage::Vector< size_t> MOBILE_ATOMS
      ) const;

      //! @brief get a new stable atom type
      AtomType FindStableAtomType
      (
        const util::SiPtr< const AtomConformationalInterface> &ATOM
      ) const;

      //! @brief resolve clashes
      //! @param INFO bond length accounting info of current molecule
      //! @param MOBILE_ATOMS atoms that can be moved during perturbation
      void ResolveClashes
      (
        MolBondInfo &INFO,
        const storage::Vector< size_t> MOBILE_ATOMS
      ) const;

      //! @brief update the MolBondInfo object
      //! @param INFO bond length accounting info of current molecule
      //! @param MOBILE_ATOMS atoms that can be moved during perturbation
      void UpdateInfo
      (
        MolBondInfo &INFO,
        const storage::Vector< size_t> MOBILE_ATOMS
      ) const;

      //! @brief check if all bond lengths meet convergence criteria
      //! @param INFO bond length accounting info of the current molecule
      //! @param MOBILE_ATOMS atoms contributing bonds that can be perturbed
      //! @return true if all bond lengths are within the convergence limit
      bool CheckConvergence
      (
        const MolBondInfo &INFO,
        const storage::Vector< size_t> MOBILE_ATOMS
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class MutateBondLengths

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_MUTATE_BOND_LENGTHS_H_

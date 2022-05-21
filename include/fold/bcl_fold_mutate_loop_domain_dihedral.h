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

#ifndef BCL_FOLD_MUTATE_LOOP_DOMAIN_DIHEDRAL_H_
#define BCL_FOLD_MUTATE_LOOP_DOMAIN_DIHEDRAL_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_fold_loop_segment.h"
#include "bcl_fold_loop_segment_sequence_order.h"
#include "find/bcl_find_collector_interface.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "math/bcl_math_mutate_interface.h"
#include "storage/bcl_storage_set.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateLoopDomainDihedral
    //! @brief is for mutating the dihedral angles of a loop domain to new values.
    //! @details The residues to be mutated are collected according to a collector interface and the dihedral angles are changed
    //! according to a dihedral angle generator behind a function interface. So the methods for determining
    //! which residues to mutate and for determining their dihedral angles are flexible.
    //!
    //! @see @link example_fold_mutate_loop_domain_dihedral.cpp @endlink
    //! @author alexanns, fischea
    //! @date Sep 6, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateLoopDomainDihedral :
      public math::MutateInterface< LoopDomain>
    {

    private:

    //////////
    // data //
    //////////

      //! the method for collecting residues whose dihedral angles will be changed to new values
      util::Implementation< find::CollectorInterface< storage::List< MutationResidue>, LoopDomain> > m_ResiduesToChangeCollector;

      //! the method that determines what the dihedral angles should be set to
      util::Implementation< math::FunctionInterfaceSerializable< MutationResidue, storage::VectorND< 2, double> > > m_DihedralGenerator;

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
      MutateLoopDomainDihedral();

      //! @brief constructor taking parameters
      //! @param RESIDUES_TO_CHANGE_COLLECTOR the method for collecting residues to change their dihedral angles
      //! @param DIHEDRAL_GENERATOR the method that determines what the dihedral angles should be set to
      MutateLoopDomainDihedral
      (
        const util::ShPtr
        <
          find::CollectorInterface< storage::List< MutationResidue>, LoopDomain>
        > &RESIDUES_TO_CHANGE_COLLECTOR,
        const util::ShPtr
        <
          math::FunctionInterfaceSerializable< MutationResidue, storage::VectorND< 2, double> >
        > &DIHEDRAL_GENERATOR
      );

      //! @brief Clone function
      //! @return pointer to new MutateLoopDomainDihedral
      MutateLoopDomainDihedral *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator taking an ARGUMENT and returning a mutated object of t_ArgumentType
      //! @param LOOP_DOMAIN Argument of interest
      //! @return MutateResult that results from mutating to the argument
      math::MutateResult< LoopDomain> operator()( const LoopDomain &LOOP_DOMAIN) const;

    //////////////////////
    // input and output //
    //////////////////////

    }; // class MutateLoopDomainDihedral

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_LOOP_DOMAIN_DIHEDRAL_H_ 

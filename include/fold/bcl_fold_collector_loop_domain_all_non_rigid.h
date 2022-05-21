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

#ifndef BCL_FOLD_COLLECTOR_LOOP_DOMAIN_ALL_NON_RIGID_H_
#define BCL_FOLD_COLLECTOR_LOOP_DOMAIN_ALL_NON_RIGID_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "find/bcl_find_collector_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CollectorLoopDomainAllNonRigid
    //! @brief is for collecting all residue in a loop domain which are not specified as rigid.
    //!
    //! @see @link example_fold_collector_loop_domain_all_non_rigid.cpp @endlink
    //! @author alexanns
    //! @date Sep 5, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CollectorLoopDomainAllNonRigid :
      public find::CollectorInterface< storage::List< MutationResidue>, LoopDomain>
    {

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
      CollectorLoopDomainAllNonRigid();

      //! @brief Clone function
      //! @return pointer to new MutateLoopDomainAllNonRigid
      CollectorLoopDomainAllNonRigid *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const
      {
        static const std::string s_alias( "CollectorLoopDomainAllNonRigid");
        return s_alias;
      }

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        io::Serializer serializer;
        serializer.SetClassDescription( "Collects non-rigid loop domains");

        return serializer;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief Collect taking an ARGUMENT and returning a list of MutationResidues
      //! @param LOOP_DOMAIN LoopDomain of interest
      //! @return storage::List< MutationResidue> that results from collecting from "LOOP_DOMAIN"
      storage::List< MutationResidue> Collect( const LoopDomain &LOOP_DOMAIN) const;

    //////////////////////
    // input and output //
    //////////////////////

    }; // class CollectorLoopDomainAllNonRigid

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_COLLECTOR_LOOP_DOMAIN_ALL_NON_RIGID_H_

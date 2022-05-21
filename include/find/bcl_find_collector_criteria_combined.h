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

#ifndef BCL_FIND_COLLECTOR_CRITERIA_COMBINED_H_
#define BCL_FIND_COLLECTOR_CRITERIA_COMBINED_H_

// include the namespace header
#include "bcl_find.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_find_collector_criteria_interface.h"
#include "assemble/bcl_assemble_sse.h"
#include "io/bcl_io_serialization.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_binary_function_interface_serializable.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace find
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CollectorCriteriaCombined
    //! @brief This class allows to collect sets of arguments based on a number of criteria.
    //! @details This class allows to pick a subset of arguments from a set of arguments that all fulfill a set of
    //!  specified criteria.
    //!
    //! @see @link example_find_collector_criteria_combined.cpp @endlink
    //! @author linders, karakam
    //! @date Aug 13, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType>
    class CollectorCriteriaCombined :
      public CollectorCriteriaInterface< util::SiPtrList< const t_ArgumentType>, util::SiPtrList< const t_ArgumentType>, t_ArgumentType>
    {
    private:

    //////////
    // data //
    //////////

      //! list of criteria used in the collect function
      //! every criterion is a BinaryFunctionInterfaceSerializable that takes two arguments and returns bool
      storage::Vector< util::Implementation< util::BinaryFunctionInterfaceSerializable< t_ArgumentType, t_ArgumentType, bool> > > m_CombinedCriteria;

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
      CollectorCriteriaCombined() :
        m_CombinedCriteria()
      {
      }

      //! @brief constructor from a single criteria
      //! @param SP_CRITERIA criteria
      CollectorCriteriaCombined
      (
        const util::ShPtr< util::BinaryFunctionInterfaceSerializable< t_ArgumentType, t_ArgumentType, bool> > &SP_CRITERIA
      ) :
        m_CombinedCriteria( 1, *SP_CRITERIA)
      {
      }

      //! @brief constructor from a list of criteria
      //! @param COMBINED_CRITERIA list of criteria used in the collect function
      CollectorCriteriaCombined
      (
        const util::ShPtrList< util::BinaryFunctionInterfaceSerializable< t_ArgumentType, t_ArgumentType, bool> > &COMBINED_CRITERIA
      ) :
        m_CombinedCriteria()
      {
        for( auto crit_it( COMBINED_CRITERIA.Begin()); crit_it != COMBINED_CRITERIA.End(); ++crit_it)
        {
          m_CombinedCriteria.PushBack( **crit_it);
        }
      }

      //! @brief Clone function
      //! @return pointer to new CollectorSSEsCriteriaCombined
      CollectorCriteriaCombined< t_ArgumentType> *Clone() const
      {
        return new CollectorCriteriaCombined< t_ArgumentType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }
      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const
      {
        static const std::string s_alias( "CollectorCriteriaCombined");
        return s_alias;
      }

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        io::Serializer serializer;
        serializer.SetClassDescription( "Combine multiple criteria.");
        serializer.AddInitializer
        (
          "combined criteria",
          "combined criteria",
          io::Serialization::GetAgent( &m_CombinedCriteria)
        );

        return serializer;
      }

    ////////////////
    // operations //
    ////////////////

      //! Collect the sses from SSES that fulfill certain criteria (see member variable) with respect to CRITERIA_SSE
      //! @param ARGUMENT_LIST List of arguments
      //! @param CRITERIA_ARGUMENT argument that will be used to initialize criteria
      //! @return returns the subset of sses from arguments that meet criteria with respect to CRITERIA_ARGUMENT
      util::SiPtrList< const t_ArgumentType> Collect
      (
        const util::SiPtrList< const t_ArgumentType> &ARGUMENT_LIST,
        const t_ArgumentType &CRITERIA_ARGUMENT
      ) const
      {
        // initialize empty SiPtrList to be returned
        util::SiPtrList< const t_ArgumentType> return_arguments;

        // iterate over the sses in SSES
        for
        (
          auto arg_itr( ARGUMENT_LIST.Begin()), arg_itr_end( ARGUMENT_LIST.End());
          arg_itr != arg_itr_end;
          ++arg_itr
        )
        {
          // initialize bool that keeps track of whether this argument satisfies all criteria
          bool all_criteria_fulfilled( true);

          // iterate over the criteria
          for
          (
            auto criteria_itr( m_CombinedCriteria.Begin()), criteria_itr_end( m_CombinedCriteria.End());
            criteria_itr != criteria_itr_end;
            ++criteria_itr
          )
          {
            // if this criteria is not satisfied
            if( !( *criteria_itr)->operator ()( **arg_itr, CRITERIA_ARGUMENT))
            {
              // set the boolean and break out of the loop
              all_criteria_fulfilled = false;
              continue;
            }
          }

          // if all criteria have been satisfied insert this argument into return list
          if( all_criteria_fulfilled)
          {
            return_arguments.PushBack( **arg_itr);
          }
        }

        // end
        return return_arguments;
      } // Collect

    //////////////////////
    // input and output //
    //////////////////////

    }; // class CollectorCriteriaCombined

    // instantiate s_Instance
    template< typename t_ArgumentType>
    const util::SiPtr< const util::ObjectInterface> CollectorCriteriaCombined< t_ArgumentType>::s_Instance
    (
      util::Enumerated< CollectorCriteriaInterface< util::SiPtrList< const t_ArgumentType>, util::SiPtrList< const t_ArgumentType>, t_ArgumentType> >::AddInstance( new CollectorCriteriaCombined< t_ArgumentType>())
    );

  } // namespace find
} // namespace bcl

#endif // BCL_FIND_COLLECTOR_CRITERIA_COMBINED_H_

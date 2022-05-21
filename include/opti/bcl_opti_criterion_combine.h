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

#ifndef BCL_OPTI_CRITERION_COMBINE_H_
#define BCL_OPTI_CRITERION_COMBINE_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opti_criterion_interface.h"
#include "storage/bcl_storage_list.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CriterionCombine
    //! @brief is a collection of CriterionInterface derived classes
    //! @details CriterionCombine stores one or more CriterionInterface derived classes and provides a single function
    //! call to let evaluate if any of the termination criteria are met.
    //!
    //! @tparam t_ArgumentType is the type of the approximation argument
    //! @tparam t_ResultType is the type of the approximation result
    //!
    //! @see @link example_opti_criterion_combine.cpp @endlink
    //! @author fischea
    //! @date Dec 13, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class CriterionCombine :
      public CriterionInterface< t_ArgumentType, t_ResultType>
    {

    //////////
    // data //
    //////////

    private:

      //! vector of ShPtrs containing the criteria combined
      storage::List< util::Implementation< CriterionInterface< t_ArgumentType, t_ResultType> > > m_CriteriaList;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      CriterionCombine() :
        m_CriteriaList()
      {
      }

      //! @brief construct from criteria list
      //! @param CRITERIA list of the termination criteria to be combined
      CriterionCombine
      (
        const storage::List< util::Implementation< CriterionInterface< t_ArgumentType, t_ResultType> > > &CRITERIA_LIST
      ) :
        m_CriteriaList( CRITERIA_LIST)
      {
      }

      //! @brief Clone function
      //! @return pointer to a new CriterionCombine< t_ArgumentType, t_ResultType>
      CriterionCombine *Clone() const
      {
        return new CriterionCombine< t_ArgumentType, t_ResultType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief returns the class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetAlias() const
      {
        static const std::string s_alias( "Any");
        return s_alias;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief adds a criterion to the list of criteria to be combined
      //! @param CRITERION criterion to be added
      void InsertCriteria( const CriterionInterface< t_ArgumentType, t_ResultType> &CRITERION)
      {
        m_CriteriaList.PushBack( CRITERION);
      }

      //! @brief returns if any of the combined termination criteria are met
      //! @param TRACKER tracker to evaluate the criterion for
      //! @return true, if any of the termination criteria are met
      bool CriteriaMet( const Tracker< t_ArgumentType, t_ResultType> &TRACKER) const
      {
        // iterate over the criteria list
        for
        (
          typename storage::List< util::Implementation< CriterionInterface< t_ArgumentType, t_ResultType> > >::const_iterator
            itr( m_CriteriaList.Begin()), itr_end( m_CriteriaList.End());
          itr != itr_end;
          ++itr
        )
        {
          if( ( *itr)->CriteriaMet( TRACKER))
          {
            return true;
          }
        }

        // only reached if none of the criteria is met
        return false;
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // read members
        io::Serialize::Read( m_CriteriaList, ISTREAM);

        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write members
        io::Serialize::Write( m_CriteriaList, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        io::Serializer serializer;
        serializer.SetClassDescription( "Triggers when any of the internal criteria are met");
        serializer.AddInitializer
        (
          "",
          "",
          io::Serialization::GetAgent( &m_CriteriaList)
        );
        return serializer;
      }

    }; // template class CriterionCombine< t_ArgumentType, t_ResultType>

    // instantiate s_Instance
    template< typename t_ArgumentType, typename t_ResultType>
    const util::SiPtr< const util::ObjectInterface> CriterionCombine< t_ArgumentType, t_ResultType>::s_Instance
    (
      util::Enumerated< CriterionInterface< t_ArgumentType, t_ResultType> >::AddInstance
      (
        new CriterionCombine< t_ArgumentType, t_ResultType>()
      )
    );

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_CRITERION_COMBINE_H_

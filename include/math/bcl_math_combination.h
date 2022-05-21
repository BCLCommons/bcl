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

#ifndef BCL_MATH_COMBINATION_H_
#define BCL_MATH_COMBINATION_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_set.h"
#include "storage/bcl_storage_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Combination
    //! @brief class for representing combinations
    //! @details handles combinations of a given size from a set
    //!
    //! @tparam t_KeyType key type used by data stored in the set
    //! @tparam t_KeyCompare compare type used by data stored in the set
    //!
    //! @see @link example_math_combination.cpp @endlink
    //! @author weinerbe
    //! @date Sep 16, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_KeyType, typename t_KeyCompare>
    class Combination :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! set containing the data
      storage::Set< t_KeyType, t_KeyCompare> m_Data;

      //! size of the combination
      size_t m_CombinationSize;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Combination() :
        m_Data(),
        m_CombinationSize( util::GetUndefinedSize_t())
      {
      }

      //! @brief construct from the set and combination sizes
      //! @param DATA set containing the elements
      //! @param COMBINATION_SIZE size of the combination
      Combination
      (
        const storage::Set< t_KeyType, t_KeyCompare> &DATA,
        const size_t COMBINATION_SIZE
      ) :
        m_Data( DATA),
        m_CombinationSize( COMBINATION_SIZE)
      {
      }

      //! @brief Clone function
      //! @return pointer to new Combination
      Combination< t_KeyType, t_KeyCompare> *Clone() const
      {
        return new Combination< t_KeyType, t_KeyCompare>( *this);
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

      //! @brief gets the data
      //! @return the data
      const storage::Set< t_KeyType, t_KeyCompare> &GetData() const
      {
        return m_Data;
      }

      //! @brief gets the combination size
      //! @return the combination size
      const size_t &GetCombinationSize() const
      {
        return m_CombinationSize;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief calculate the total possible number of combinations
      //! @return the total possible number of combinations
      size_t GetNumberOfCombinations() const
      {
        // assert that the combination size is smaller than the set size
        BCL_Assert
        (
          m_CombinationSize <= m_Data.GetSize(),
          "The size of the combination cannot exceed the size of the set"
        );

        // calculate the number of combinations via the binomial coefficient
        return BinomialCoefficient( m_Data.GetSize(), m_CombinationSize);
      }

      //! @brief returns a list of all possible combinations
      //! @return list of all possible combinations
      storage::List< storage::Set< t_KeyType, t_KeyCompare> > GetAllCombinations() const
      {
        // assert that the combination size is smaller than the set size
        BCL_Assert
        (
          m_CombinationSize <= m_Data.GetSize(),
          "The size of the combination cannot exceed the size of the set"
        );

        // initialize list of combinations
        storage::List< storage::Set< t_KeyType, t_KeyCompare> > combinations;

        // initialize position vector
        storage::Vector< size_t> positions( m_CombinationSize);

        // get the combinations
        GetCombinations( positions, 0, 0, combinations);

        // end
        return combinations;
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
        io::Serialize::Read( m_Data           , ISTREAM);
        io::Serialize::Read( m_CombinationSize, ISTREAM);

        // return the stream
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write members
        io::Serialize::Write( m_Data           , OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_CombinationSize, OSTREAM, INDENT);

        // return the stream
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief recursively calculates all combinations
      //! @param POSITIONS vector containing current indices to build the combination
      //! @param DEPTH how many elements have already been chosen
      //! @param MARGIN the first position to look for a new element (left margin)
      //! @param COMBINATION_LIST list containing all calculated combinations
      void GetCombinations
      (
        storage::Vector< size_t> &POSITIONS,
        const size_t &DEPTH,
        const size_t &MARGIN,
        storage::List< storage::Set< t_KeyType, t_KeyCompare> > &COMBINATION_LIST
      ) const
      {
        // if the combination is the right size
        if( DEPTH >= m_CombinationSize)
        {
          // create a set to store the combination
          storage::Set< t_KeyType, t_KeyCompare> combination;

          // iterate through the positions
          for
          (
            storage::Vector< size_t>::const_iterator itr( POSITIONS.Begin()), itr_end( POSITIONS.End());
            itr != itr_end;
            ++itr
          )
          {
            // get an iterator on the data
            typename storage::Set< t_KeyType, t_KeyCompare>::const_iterator set_itr( m_Data.Begin());

            // move to the position given by "itr"
            for( size_t i( 0), i_end( *itr); i != i_end; ++i)
            {
              ++set_itr;
            }

            // pushback the geometry at that position into the combination
            combination.InsertElement( *set_itr);
          }

          // add the combination to the vector
          COMBINATION_LIST.PushBack( combination);

          // end
          return;
        }

        // exit if there are not enough remaining elements
        if( m_Data.GetSize() - MARGIN < m_CombinationSize - DEPTH)
        {
          return;
        }

        // select new elements to the right of the previous selection
        for( size_t i( MARGIN); i < m_Data.GetSize(); ++i)
        {
          POSITIONS( DEPTH) = i;
          GetCombinations( POSITIONS, DEPTH + 1, i + 1, COMBINATION_LIST);
        }
      }

    }; // template class Combination

    // instantiate s_Instance
    template< typename t_KeyType, typename t_KeyCompare>
    const util::SiPtr< const util::ObjectInterface> Combination< t_KeyType, t_KeyCompare>::s_Instance
    (
      GetObjectInstances().AddInstance( new Combination< t_KeyType, t_KeyCompare>())
    );

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_COMBINATION_H_

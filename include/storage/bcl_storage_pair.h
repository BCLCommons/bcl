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

#ifndef BCL_STORAGE_PAIR_H_
#define BCL_STORAGE_PAIR_H_

// include the namespace header
#include "bcl_storage.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialize.h"
#include "util/bcl_util_binary_function_interface.h"
#include "util/bcl_util_object_instances.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace storage
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Pair
    //! @brief A helper class, analogous to std::pair.
    //! @details It allows to combine two different data types without writing an own class or class for it.
    //!
    //! attributes:
    //! - holds two objects of possibly different data type
    //! - no iterators
    //! - no insertion or removal
    //! - direct access to either of the two objects is available through the First() and Second() functions
    //!
    //! @tparam t_First the type of the first object
    //! @tparam t_Second the type of the second object
    //!
    //! @see @link example_storage_pair.cpp @endlink
    //! @author woetzen, haenigc
    //! @date 1/2/2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_First, typename t_Second>
    class Pair :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      t_First m_First;   //!< First object
      t_Second m_Second; //!< Second object

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
      Pair() :
        m_First(),
        m_Second()
      {
      }

      //! @brief construct Pair from two values
      //! @param FIRST the object which will be m_First
      //! @param SECOND the object which will be m_Second
      Pair( const t_First &FIRST, const t_Second &SECOND) :
        m_First( FIRST),
        m_Second( SECOND)
      {
      }

      //! @brief construct Pair from std::pair
      //! @param PAIR which will be the Pair
      Pair( const std::pair< t_First, t_Second> &PAIR) :
        m_First( PAIR.first),
        m_Second( PAIR.second)
      {
      };

      //! copy constructor
      Pair( const Pair &PAIR) :
        m_First( PAIR.m_First),
        m_Second( PAIR.m_Second)
      {
      }

      //! move constructor
      Pair( Pair && PAIR) :
        m_First( std::move( PAIR.m_First)),
        m_Second( std::move( PAIR.m_Second))
      {
      }

      //! virtual copy constructor (Clone)
      Pair< t_First, t_Second> *Clone() const
      {
        return new Pair< t_First, t_Second>( *this);
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

      //! @brief First gives access to the first object
      //! @return returns m_First
      t_First &First()
      {
        return m_First;
      }

      //! @brief First gives const access to the first object
      //! @return returns const reference to m_First
      t_First const &First() const
      {
        return m_First;
      }

      //! @brief Second gives access to the second object
      //! @return returns reference to m_Second
      t_Second &Second()
      {
        return m_Second;
      }

      //! @brief Second gives const access to the second object
      //! @return returns const reference to m_Second
      t_Second const &Second() const
      {
        return m_Second;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator = defines equating two Pairs
      //! @param PAIR which the current Pair will be set equal to
      //! @return returns Pair which as been set equal to PAIR
      Pair &operator =( const Pair &PAIR)
      {
        m_First  = PAIR.m_First;  //< set first object equal to First object of PAIR
        m_Second = PAIR.m_Second; //< set second object equal to Second object of Pair

        // return the equated Pair
        return *this;
      }

      //! @brief operator = defines moving two Pairs
      //! @param PAIR which the current Pair will be set equal to
      //! @return returns Pair which as been set equal to PAIR
      Pair &operator =( Pair && PAIR)
      {
        if( this != &PAIR)
        {
          m_First  = std::move( PAIR.m_First);  //< set first object equal to First object of PAIR
          m_Second = std::move( PAIR.m_Second); //< set second object equal to Second object of Pair
        }

        // return the equated Pair
        return *this;
      }

      //! @brief operator = defines equating a Pair to a std::pair
      //! @param PAIR which the current Pair will be set equal to
      //! @return returns Pair which as been set equal to PAIR
      Pair< t_First, t_Second> &operator =( const std::pair< t_First, t_Second> &PAIR)
      {
        m_First = PAIR.first; //< set first object equal to First object of PAIR
        m_Second = PAIR.second; //< set second object equal to Second object of PAIR
        return *this; //< return the equated Pair
      }

      //! @brief operator < for checking if one pair is less than the other one
      //! true if first of PAIR_A is less than first of PAIR_B or if the first of both pairs
      //! is the same and the second of PAIR_A is less than the second of PAIR_B
      //! @param PAIR_A first pair for comparison
      //! @param PAIR_B second pair for comparison
      //! @return comparison result of PAIR_A and PAIR_B
      bool operator<( const Pair< t_First, t_Second> &PAIR_B) const
      {
        return m_First < PAIR_B.m_First || ( m_First == PAIR_B.m_First && m_Second < PAIR_B.m_Second);
      }

      //! @brief operator > for checking if one pair is less than the other one
      //! true if first of PAIR_A is less than first of PAIR_B or if the first of both pairs
      //! is the same and the second of PAIR_A is less than the second of PAIR_B
      //! @param PAIR_A first pair for comparison
      //! @param PAIR_B second pair for comparison
      //! @return comparison result of PAIR_A and PAIR_B
      bool operator>( const Pair< t_First, t_Second> &PAIR_B) const
      {
        return m_First > PAIR_B.m_First || ( m_First == PAIR_B.m_First && m_Second > PAIR_B.m_Second);
      }

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! @brief Read inputs a Pair from a stream
      //! @param ISTREAM is the input stream
      //! @return returns the input stream
      std::istream &Read( std::istream &ISTREAM)
      {
        // read members
        io::Serialize::Read( m_First, ISTREAM);
        io::Serialize::Read( m_Second, ISTREAM);

        // end
        return ISTREAM;
      }

      //! @brief Write outputs Pair to a stream
      //! @param OSTREAM is the output stream
      //! @param INDENT indentation
      //! @return returns the output stream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write members
        io::Serialize::Write( m_First, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Second, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    }; // end Pair class

    // instantiate s_Instance
    template< typename t_First, typename t_Second>
    const util::SiPtr< const util::ObjectInterface> Pair< t_First, t_Second>::s_Instance
    (
      GetObjectInstances().AddInstance( new Pair< t_First, t_Second>())
    );

    //! @brief operator == checks if two Pairs are same
    //! @param PAIR_A the left hand operand of the operator
    //! @param PAIR_B the right hand operand of the operator
    //! @return returns a bool, true if PAIR_A and PAIR_B are the same
    template< typename t_First, typename t_Second>
    inline bool operator ==( const Pair< t_First, t_Second> &PAIR_A, const Pair< t_First, t_Second> &PAIR_B)
    {
      return ( PAIR_A.First() == PAIR_B.First() && PAIR_A.Second() == PAIR_B.Second());
    }

    //! @brief operator != checks if two Pairs are different
    //! @param PAIR_A the left hand operand of the operator
    //! @param PAIR_B the right hand operand of the operator
    //! @return returns bool, true if two std::pairs are not the same
    template< typename t_First, typename t_Second>
    inline bool operator !=( const Pair< t_First, t_Second> &PAIR_A, const Pair< t_First, t_Second> &PAIR_B)
    {
      return !( PAIR_A == PAIR_B);
    }

  ////////////////////////////
  // sorting and comparison //
  ////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @class PairBinaryPredicateFirst
    //! @brief defines a binary predicate that is used to compare the first value of two pairs
    //! @remarks example unnecessary
    //! @author woetzen
    //! @date Oct 9, 2008
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_First, typename t_Second>
    class PairBinaryPredicateFirst :
      public util::BinaryFunctionInterface< Pair< t_First, t_Second>, Pair< t_First, t_Second>, bool>
    {

    //////////
    // data //
    //////////

      //! binary predicate function
      util::ShPtr< util::BinaryFunctionInterface< t_First, t_First, bool> >
        m_BinaryPredicateFirst;

    public:
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! construct from binary predicate
      PairBinaryPredicateFirst
      (
        const util::BinaryFunctionInterface< t_First, t_First, bool> &BINARY_PREDICATE
      ) :
        m_BinaryPredicateFirst( BINARY_PREDICATE.Clone())
      {
      }

      //! copy constructor
      PairBinaryPredicateFirst< t_First, t_Second> *Clone() const
      {
        return new PairBinaryPredicateFirst< t_First, t_Second>( *this);
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

    //////////////
    // operator //
    //////////////

      //! @brief operator () for comparing two pairs by their first elements
      //! @param PAIR_A first pair for comparison
      //! @param PAIR_B second pair for comparison
      //! @return boolean value, true if two pairs has the same first element
      bool operator()( const Pair< t_First, t_Second> &PAIR_A, const Pair< t_First, t_Second> &PAIR_B) const
      {
        // compare first elements and return result
        return m_BinaryPredicateFirst->operator()( PAIR_A.First(), PAIR_B.First());
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
        // read member
        io::Serialize::Read( m_BinaryPredicateFirst, ISTREAM);

        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write member
        io::Serialize::Write( m_BinaryPredicateFirst, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    }; // template class PairBinaryPredicateFirst

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @class PairBinaryPredicateSecond
    //! @brief defines a binary predicate that is used to compare the second value of two pairs
    //! @remarks example unnecessary
    //! @author woetzen
    //! @date Oct 9, 2008
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_First, typename t_Second>
    class PairBinaryPredicateSecond :
      public util::BinaryFunctionInterface< Pair< t_First, t_Second>, Pair< t_First, t_Second>, bool>
    {
    //////////
    // data //
    //////////

      //! binary predicate function
      util::ShPtr< util::BinaryFunctionInterface< t_Second, t_Second, bool> >
        m_BinaryPredicateSecond;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! construct from binary predicate
      PairBinaryPredicateSecond
      (
        const util::BinaryFunctionInterface< t_Second, t_Second, bool> &BINARY_PREDICATE
      ) :
        m_BinaryPredicateSecond( BINARY_PREDICATE.Clone())
      {
      }

      //! copy constructor
      PairBinaryPredicateSecond< t_First, t_Second> *Clone() const
      {
        return new PairBinaryPredicateSecond< t_First, t_Second>( *this);
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

    //////////////
    // operator //
    //////////////

      //! @brief operator () for comparing two pairs by their second elements
      //! @param PAIR_A first pair for comparison
      //! @param PAIR_B second pair for comparison
      //! @return boolean value, true if two pairs has the same first element
      bool operator()( const Pair< t_First, t_Second> &PAIR_A, const Pair< t_First, t_Second> &PAIR_B) const
      {
        // compare first elements and return result
        return m_BinaryPredicateSecond->operator()( PAIR_A.Second(), PAIR_B.Second());
      }

      //! @brief operator () for comparing two pairs by their second elements
      //! @param PAIR_A first pair for comparison
      //! @param PAIR_B second pair for comparison
      //! @return boolean value, true if two pairs has the same first element
      bool operator()( const std::pair< t_First, t_Second> &PAIR_A, const std::pair< t_First, t_Second> &PAIR_B) const
      {
        // compare first elements and return result
        return m_BinaryPredicateSecond->operator()( PAIR_A.second(), PAIR_B.second());
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
        // read member
        io::Serialize::Read( m_BinaryPredicateSecond, ISTREAM);

        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write member
        io::Serialize::Write( m_BinaryPredicateSecond, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    }; // template class PairBinaryPredicateSecond

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @class PairEqualFirst
    //! @brief checks if a pair has the given first element
    //! @details Needed to order use std::find_if and similar functions for pairs with given first element
    //! @remarks example unnecessary
    //! @author woetzen
    //! @date Oct 9, 2008
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_First>
    class PairEqualFirst
    {
      //! first element to compare
      const t_First m_First;

    public:

      //! constructor from a type t_First
      PairEqualFirst( const t_First &FIRST) :
        m_First( FIRST)
      {
      }

      //! @brief operator () for checking if a pair has the given first element
      //! @param PAIR pair for comparison
      //! @return boolean value, true if two pairs has the same first element
      template< typename t_Second>
      bool operator()( const Pair< t_First, t_Second> &PAIR) const
      {
        // compare first elements and return
        return PAIR.First() == m_First;
      }

    };

  } // namespace storage
} // namespace bcl

#endif // BCL_STORAGE_PAIR_H_

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

#ifndef BCL_STORAGE_TRIPLET_H_
#define BCL_STORAGE_TRIPLET_H_

// include the namespace header
#include "bcl_storage.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialize.h"
#include "util/bcl_util_object_instances.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace storage
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Triplet
    //! @brief A helper class, in analogy to std::pair.
    //! @details It allows to combine three different data types without writing an own class or struct for it.
    //!
    //! attributes:
    //! - holds three objects of possibly different data types
    //! - no iterators
    //! - no insertion or removal
    //! - direct access to either of the three objects is available through the First() and Second() and Third()
    //!
    //! @tparam t_First the type of the first object
    //! @tparam t_Second the type of the second object
    //! @tparam t_Third the type of the third object
    //!
    //! @see @link example_storage_triplet.cpp @endlink
    //! @author staritrd, meilerj
    //! @date 2/1/2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_First, typename t_Second, typename t_Third>
    class Triplet :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      t_First m_First;    //!< First object
      t_Second m_Second;   //!< Second object
      t_Third m_Third;    //!< Third object

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Triplet() :
        m_First(), m_Second(), m_Third()
      {
      }

      //! @brief construct Triplet from three values
      //! @param FIRST the first object to go into the triplet
      //! @param SECOND the second object to go into the triplet
      //! @param THIRD the third object to go into the triplet
      Triplet( const t_First &FIRST, const t_Second &SECOND, const t_Third &THIRD) :
        m_First( FIRST),
        m_Second( SECOND),
        m_Third( THIRD)
      {
      }

      //! @brief copy constructor
      //! @param TRIPLET Triplet for which a copy is desired
      Triplet( const Triplet &TRIPLET) :
        m_First( TRIPLET.m_First),
        m_Second( TRIPLET.m_Second),
        m_Third( TRIPLET.m_Third)
      {
      }

      //! @brief virtual copy constructor, needed to make hard copies from pointer to it
      virtual Triplet *Clone() const
      {
        return new Triplet( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      virtual const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief First gives the first value in the triplet
      //! @return gives a changable reference to the first member
      virtual t_First &First()
      {
        return m_First;
      }

      //! @brief First gives the first value in the triplet as a constant reference
      //! @return gives the first value in the triplet as a constant reference
      virtual t_First const &First() const
      {
        return m_First;
      }

      //! @brief Second gives the second value in the triplet
      //! @brief gives the second value in the triplet as a changable reference
      virtual t_Second &Second()
      {
        return m_Second;
      }

      //! @brief Second gives the second value in the triplet as a constant reference
      //! @return gives the second value in the triplet as a constant reference
      virtual t_Second const &Second() const
      {
        return m_Second;
      }

      //! @brief Third gives the third value in the triplet
      //! @brief gives the third value in the triplet as a changable reference
      virtual t_Third &Third()
      {
        return m_Third;
      }

      //! @brief Third gives the third value in the triplet as a constant reference
      //! @return gives the third value in the triplet as a constant reference
      virtual t_Third const &Third() const
      {
        return m_Third;
      }

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! @brief Write outputs Triplet to a stream
      //! @param OSTREAM is the output stream
      //! @param INDENT indentation
      //! @return returns the output stream
      virtual std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write members
        io::Serialize::Write( m_First, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Second, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Third, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

      //! @brief Read inputs a Triplet from a stream
      //! @param ISTREAM is the input stream
      //! @return returns the input stream
      virtual std::istream &Read( std::istream &ISTREAM)
      {
        // read members
        io::Serialize::Read( m_First, ISTREAM);
        io::Serialize::Read( m_Second, ISTREAM);
        io::Serialize::Read( m_Third, ISTREAM);

        // end
        return ISTREAM;
      }

    }; // end Triplet class

    // instantiate s_Instance
    template< typename t_First, typename t_Second, typename t_Third>
    const util::SiPtr< const util::ObjectInterface> Triplet< t_First, t_Second, t_Third>::s_Instance
    (
      GetObjectInstances().AddInstance( new Triplet< t_First, t_Second, t_Third>())
    );

    //! @brief operator < checks which of two Triplets is the "larger"
    //! @param TRIPLET_A left hand operand and first Triplet for comparison
    //! @param TRIPLET_B right hand operand and second Triplet for comparison
    //! Needed to order std::set (e.g. std::set< Triplet< t_First, t_Second, t_Third> >)
    template< typename t_First, typename t_Second, typename t_Third>
    inline bool operator <
    (
      const Triplet< t_First, t_Second, t_Third> &TRIPLET_A,
      const Triplet< t_First, t_Second, t_Third> &TRIPLET_B
    )
    {
      // return true if first member of TRIPLET_A is less than first member of TRIPLET_B
      if( TRIPLET_A.First() < TRIPLET_B.First())
      {
        return true;
      }
      if( TRIPLET_A.First() > TRIPLET_B.First())
      {
        return false;
      }

      // return true if second member of TRIPLET_A is less than second member of TRIPLET_B
      if( TRIPLET_A.Second() < TRIPLET_B.Second())
      {
        return true;
      }
      if( TRIPLET_B.Second() < TRIPLET_A.Second())
      {
        return false;
      }

      // return true if third member of TRIPLET_A is less than third member of TRIPLET_B
      if( TRIPLET_A.Third() < TRIPLET_B.Third())
      {
        return true;
      }

      // if no member of TRIPLET_A is less than the corresponding member of TRIPLET_B then return false
      return false;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LessThanFirst
    //! @brief binary operator function struct that only compares the first element of the triplet
    //! @author mendenjl
    //! @date Apr 09, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API LessThanFirst
    {
    public:

      template< typename t_First, typename t_Second, typename t_Third>
      bool operator()
      (
        const Triplet< t_First, t_Second, t_Third> &TRIPLET_A,
        const Triplet< t_First, t_Second, t_Third> &TRIPLET_B
      ) const
      {
        // compare first element only
        return TRIPLET_A.First() < TRIPLET_B.First();
      }
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LessThanSecond
    //! @brief binary operator function struct that only compares the second element of the triplet
    //! @author mendenjl
    //! @date Apr 09, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API LessThanSecond
    {
    public:

      template< typename t_First, typename t_Second, typename t_Third>
      bool operator()
      (
        const Triplet< t_First, t_Second, t_Third> &TRIPLET_A,
        const Triplet< t_First, t_Second, t_Third> &TRIPLET_B
      ) const
      {
        // compare 2nd element only
        return TRIPLET_A.Second() < TRIPLET_B.Second();
      }
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LessThanThird
    //! @brief binary operator function struct that only compares the third element of the triplet
    //! @author mendenjl
    //! @date Apr 09, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API LessThanThird
    {
    public:

      template< typename t_First, typename t_Second, typename t_Third>
      bool operator()
      (
        const Triplet< t_First, t_Second, t_Third> &TRIPLET_A,
        const Triplet< t_First, t_Second, t_Third> &TRIPLET_B
      ) const
      {
        // compare 3rd element only
        return TRIPLET_A.Third() < TRIPLET_B.Third();
      }
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class GreaterThanThird
    //! @brief binary operator function struct that only compares the third element of the triplet
    //! @author mendenjl
    //! @date Apr 09, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API GreaterThanThird
    {
    public:

      template< typename t_First, typename t_Second, typename t_Third>
      bool operator()
      (
        const Triplet< t_First, t_Second, t_Third> &TRIPLET_A,
        const Triplet< t_First, t_Second, t_Third> &TRIPLET_B
      ) const
      {
        // compare 3rd element only
        return TRIPLET_A.Third() > TRIPLET_B.Third();
      }
    };

    //! @brief operator == tests wether two Triplets are the same
    //! @param TRIPLET_A left hand operand and first Triplet for comparison
    //! @param TRIPLET_B right hand operand and second Triplet for comparison
    template< typename t_First, typename t_Second, typename t_Third>
    inline bool operator ==
    (
      const Triplet< t_First, t_Second, t_Third> &TRIPLET_A,
      const Triplet< t_First, t_Second, t_Third> &TRIPLET_B
    )
    {
      // return true if all members of both Triplets are identical
      return
        TRIPLET_A.First() == TRIPLET_B.First() &&
        TRIPLET_A.Second() == TRIPLET_B.Second() &&
        TRIPLET_A.Third() == TRIPLET_B.Third();
    }

  } // namespace storage
} // namespace bcl

#endif //BCL_STORAGE_TRIPLET_H_

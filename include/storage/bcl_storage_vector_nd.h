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

#ifndef BCL_STORAGE_VECTOR_ND_H_
#define BCL_STORAGE_VECTOR_ND_H_

// include the namespace header
#include "bcl_storage.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_storage_vector.h"
#include "util/bcl_util_object_instances.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace storage
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class VectorND
    //! @brief template class VectorND< N, t_DataType>
    //! @details attributes:
    //! a.) has a constant number of elements it holds
    //! b.) no iterators exist for this container
    //! c.) random access to elements
    //! d.) no insertion or deletion of elements
    //!
    //! @see @link example_storage_vector_nd.cpp @endlink
    //! @author woetzen
    //! @date Oct 29, 2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< unsigned int N, typename t_DataType>
    class VectorND :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      Vector< t_DataType> m_Data; //! this is the complete data

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      VectorND() :
        m_Data( N)
      {
      }

      //! virtual copy constructor
      VectorND< N, t_DataType> *Clone() const
      {
        return new VectorND< N, t_DataType>( *this);
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

      //! return const reference to first element
      t_DataType const &First() const
      {
        return m_Data( 0);
      }

      //! return reference to changeable element
      t_DataType       &First()
      {
        return m_Data( 0);
      }

      //! return const reference to second element
      t_DataType const &Second() const
      {
        return m_Data( 1);
      }

      //! return reference to second element
      t_DataType       &Second()
      {
        return m_Data( 1);
      }

      //! return const reference to third element
      t_DataType const &Third() const
      {
        return m_Data( 2);
      }

      //! return reference to third element
      t_DataType       &Third()
      {
        return m_Data( 2);
      }

      //! return const reference to fourth element
      t_DataType const &Fourth() const
      {
        return m_Data( 3);
      }

      //! return reference to fourth element
      t_DataType       &Fourth()
      {
        return m_Data( 3);
      }

    ///////////////
    // operators //
    ///////////////

      //! operator( POS) return const reference to element at POS
      t_DataType const &operator()( const size_t POS) const
      {
        return m_Data( POS);
      }

      //! operator( POS) return reference to changeable element at POS
      t_DataType       &operator()( const size_t POS)
      {
        return m_Data( POS);
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! write VectorND< N, t_DataType> to std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        //write data
        io::Serialize::Write( m_Data, OSTREAM, INDENT);

        //return
        return OSTREAM;
      }

      //! read VectorND< N, t_DataType> from std::istream
      std::istream &Read( std::istream &ISTREAM)
      {
        //read data
        io::Serialize::Read( m_Data, ISTREAM);

        //end
        return ISTREAM;
      }

    }; // template class VectorND

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class VectorND
    //! @brief template class VectorND< 1, t_DataType>
    //! @details attributes:
    //! a.) has a constant number of elements it holds
    //! b.) no iterators exist for this container
    //! c.) random access to elements
    //! d.) no insertion or deletion of elements
    //!
    //! @see @link example_storage_vector_nd.cpp @endlink
    //! @author woetzen
    //! @date Oct 29, 2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class VectorND< 1, t_DataType> :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      t_DataType m_First; //! this is the data for the first element;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      VectorND()
      {
      }

      //! construct from one element FIRST
      VectorND( const t_DataType &FIRST) :
        m_First( FIRST)
      {
      }

      //! virtual copy constructor
      VectorND< 1, t_DataType> *Clone() const
      {
        return new VectorND< 1, t_DataType>( *this);
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

      //! return const refernce to first element
      t_DataType const &First() const
      {
        return m_First;
      }

      //! return reference to changable element
      t_DataType       &First()
      {
        return m_First;
      }

    ///////////////
    // operators //
    ///////////////

      //! operator( POS) return const refernce to element at POS
      t_DataType const &operator()( const size_t POS) const
      {
        switch( POS)
        {
          //first element
          case 0:
          {
            return m_First;
          }

          //if pos is higher than number elements
          default:
          {
            BCL_Exit( "try to access non existing element in VectorND< 1, t_DataType>! " + util::Format()( POS) + " >= 1", -1);
            return m_First;
          }
        }
      }

      //! operator( POS) return reference to changable element at POS
      t_DataType       &operator()( const size_t POS)
      {
        switch( POS)
        {
          //first element
          case 0:
          {
            return m_First;
          }

          //if pos is higher than number elements
          default:
          {
            BCL_Exit( "try to access non existing element in VectorND< 1, t_DataType>! " + util::Format()( POS) + " >= 1", -1);
            return m_First;
          }
        }
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! write VectorND< 1, t_DataType> to std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        //write data
        io::Serialize::Write( m_First, OSTREAM, INDENT);

        //return
        return OSTREAM;
      }

      //! read VectorND< 1, t_DataType> from std::istream
      std::istream &Read( std::istream &ISTREAM)
      {
        //read data
        io::Serialize::Read( m_First, ISTREAM);

        //end
        return ISTREAM;
      }

    }; // template class VectorND< 1, t_DataType>

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class VectorND
    //! @brief template class VectorND< 2, t_DataType>
    //! @details attributes:
    //! a.) has a constant number of elements it holds
    //! b.) no iterators exist for this container
    //! c.) random access to elements
    //! d.) no insertion or deletion of elements
    //!
    //! @see @link example_storage_vector_nd.cpp @endlink
    //! @author woetzen
    //! @date Oct 29, 2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class VectorND< 2, t_DataType> :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      t_DataType m_First;  //! this is the data for the first element;
      t_DataType m_Second; //! this is the data for the second element;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      VectorND()
      {
      }

      //! construct from one ELEMENT
      VectorND( const t_DataType &ELEMENT) :
        m_First( ELEMENT),
        m_Second( ELEMENT)
      {
      }

      //! construct from two elements FIRST and SECOND
      VectorND( const t_DataType &FIRST, const t_DataType &SECOND) :
        m_First( FIRST),
        m_Second( SECOND)
      {
      }

      //! copy constructor
      VectorND( const VectorND &A) :
        m_First( A.m_First),
        m_Second( A.m_Second)
      {
      }

      //! move constructor
      VectorND( VectorND && A) :
        m_First( std::move( A.m_First)),
        m_Second( std::move( A.m_Second))
      {
      }

      //! virtual copy constructor
      VectorND< 2, t_DataType> *Clone() const
      {
        return new VectorND< 2, t_DataType>( *this);
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

      //! return const refernce to first element
      t_DataType const &First() const
      {
        return m_First;
      }

      //! return reference to changable element
      t_DataType       &First()
      {
        return m_First;
      }

      //! return const refernce to second element
      t_DataType const &Second() const
      {
        return m_Second;
      }

      //! return reference to second element
      t_DataType       &Second()
      {
        return m_Second;
      }

    ///////////////
    // operators //
    ///////////////

      //! operator( POS) return const reference to element at POS
      t_DataType const &operator()( const size_t POS) const
      {
        switch( POS)
        {
          //first element
          case 0:
          {
            return m_First;
          }
          //second element
          case 1:
          {
            return m_Second;
          }
          //if pos is higher than number elements
          default:
          {
            BCL_Exit( "try to access non existing element in VectorND< 2, t_DataType>! " + util::Format()( POS) + " >= 2", -1);
            return m_First;
          }
        }
      }

      //! operator( POS) return reference to changeable element at POS
      t_DataType       &operator()( const size_t POS)
      {
        switch( POS)
        {
          //first element
          case 0:
          {
            return m_First;
          }
          //second element
          case 1:
          {
            return m_Second;
          }
          //if pos is higher than number elements
          default:
          {
            BCL_Exit( "try to access non existing element in VectorND< 2, t_DataType>! " + util::Format()( POS) + " >= 2", -1);
            return m_First;
          }
        }
      }

      //! @brief equal operator
      VectorND &operator =( const VectorND &VECTOR)
      {
        if( this != &VECTOR)
        {
          m_First = VECTOR.m_First;
          m_Second = VECTOR.m_Second;
        }
        return *this;
      }

      //! @brief move assignment operator
      VectorND &operator =( VectorND && VECTOR)
      {
        if( this != &VECTOR)
        {
          m_First = std::move( VECTOR.m_First);
          m_Second = std::move( VECTOR.m_Second);
        }
        return *this;
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! write VectorND< 2, t_DataType> to std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        //write data
        io::Serialize::Write( m_First, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Second, OSTREAM, INDENT);

        //return
        return OSTREAM;
      }

      //! read VectorND< 2, t_DataType> from std::istream
      std::istream &Read( std::istream &ISTREAM)
      {
        //read data
        io::Serialize::Read( m_First, ISTREAM);
        io::Serialize::Read( m_Second, ISTREAM);

        //end
        return ISTREAM;
      }

    }; // template class VectorND< 2, t_DataType>

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class VectorND
    //! @brief template class VectorND< 3, t_DataType>
    //! @details attributes:
    //! a.) has a constant number of elements it holds
    //! b.) no iterators exist for this container
    //! c.) random access to elements
    //! d.) no insertion or deletion of elements
    //!
    //! @see @link example_storage_vector_nd.cpp @endlink
    //! @author woetzen
    //! @date Oct 29, 2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class VectorND< 3, t_DataType> :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      t_DataType m_First;  //! this is the data for the first element;
      t_DataType m_Second; //! this is the data for the second element;
      t_DataType m_Third;  //! this is the data for the third element;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      VectorND()
      {
      }

      //! construct from one ELEMENT
      VectorND( const t_DataType &ELEMENT) :
        m_First( ELEMENT),
        m_Second( ELEMENT),
        m_Third( ELEMENT)
      {
      }

      //! copy constructor
      VectorND( const VectorND &A) :
        m_First( A.m_First),
        m_Second( A.m_Second),
        m_Third( A.m_Third)
      {
      }

      //! move
      VectorND( VectorND && A) :
        m_First( std::move( A.m_First)),
        m_Second( std::move( A.m_Second)),
        m_Third( std::move( A.m_Third))
      {
      }

      //! construct from three elements FIRST, SECOND and THRID
      VectorND( const t_DataType &FIRST, const t_DataType &SECOND, const t_DataType &THIRD) :
        m_First( FIRST),
        m_Second( SECOND),
        m_Third( THIRD)
      {
      }

      //! virtual copy constructor
      VectorND< 3, t_DataType> *Clone() const
      {
        return new VectorND< 3, t_DataType>( *this);
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

      //! return const refernce to first element
      t_DataType const &First() const
      {
        return m_First;
      }

      //! return reference to changable element
      t_DataType       &First()
      {
        return m_First;
      }

      //! return const refernce to second element
      t_DataType const &Second() const
      {
        return m_Second;
      }

      //! return reference to second element
      t_DataType       &Second()
      {
        return m_Second;
      }

      //! return const refernce to third element
      t_DataType const &Third() const
      {
        return m_Third;
      }

      //! return reference to third element
      t_DataType       &Third()
      {
        return m_Third;
      }

    ///////////////
    // operators //
    ///////////////

      //! operator( POS) return const reference to element at POS
      t_DataType const &operator()( const size_t POS) const
      {
        switch( POS)
        {
          //first element
          case 0:
          {
            return m_First;
          }
          //second element
          case 1:
          {
            return m_Second;
          }
          //third element
          case 2:
          {
            return m_Third;
          }
          //if pos is higher than number elements
          default:
          {
            BCL_Exit( "try to access non existing element in VectorND< 3, t_DataType>! " + util::Format()( POS) + " >= 3", -1);
            return m_First;
          }
        }
      }

      //! operator( POS) return reference to changeable element at POS
      t_DataType       &operator()( const size_t POS)
      {
        switch( POS)
        {
          //first element
          case 0:
          {
            return m_First;
          }
          //second element
          case 1:
          {
            return m_Second;
          }
          //third element
          case 2:
          {
            return m_Third;
          }
          //if pos is higher than number elements
          default:
          {
            BCL_Exit( "try to access non existing element in VectorND< 3, t_DataType>! " + util::Format()( POS) + " >= 3", -1);
            return m_First;
          }
        }
      }

      //! @brief equal operator
      VectorND &operator =( const VectorND &VECTOR)
      {
        if( this != &VECTOR)
        {
          m_First = VECTOR.m_First;
          m_Second = VECTOR.m_Second;
          m_Third = VECTOR.m_Third;
        }
        return *this;
      }

      //! @brief move assignment operator
      VectorND &operator =( VectorND && VECTOR)
      {
        if( this != &VECTOR)
        {
          m_First = std::move( VECTOR.m_First);
          m_Second = std::move( VECTOR.m_Second);
          m_Third = std::move( VECTOR.m_Third);
        }
        return *this;
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! write VectorND< 3, t_DataType> to std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        //write data
        io::Serialize::Write( m_First, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Second, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Third, OSTREAM, INDENT);

        //return
        return OSTREAM;
      }

      //! read VectorND< 3, t_DataType> from std::istream
      std::istream &Read( std::istream &ISTREAM)
      {
        //read data
        io::Serialize::Read( m_First, ISTREAM);
        io::Serialize::Read( m_Second, ISTREAM);
        io::Serialize::Read( m_Third, ISTREAM);

        //end
        return ISTREAM;
      }

    }; // template class VectorND< 3, t_DataType>

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class VectorND
    //! @brief template class VectorND< 4, t_DataType>
    //! @details attributes:
    //! a.) has a constant number of elements it holds
    //! b.) no iterators exist for this container
    //! c.) random access to elements
    //! d.) no insertion or deletion of elements
    //!
    //! @see @link example_storage_vector_nd.cpp @endlink
    //! @author woetzen
    //! @date Oct 29, 2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class VectorND< 4, t_DataType> :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      t_DataType m_First;  //! this is the data for the first element;
      t_DataType m_Second; //! this is the data for the second element;
      t_DataType m_Third;  //! this is the data for the third element;
      t_DataType m_Fourth; //! this is the data for the fourth element;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      VectorND()
      {
      }

      //! construct from one ELEMENT
      VectorND( const t_DataType &ELEMENT) :
        m_First( ELEMENT),
        m_Second( ELEMENT),
        m_Third( ELEMENT),
        m_Fourth( ELEMENT)
      {
      }

      //! construct from four elements FIRST, SECOND, THRID and FOURTH
      VectorND( const t_DataType &FIRST, const t_DataType &SECOND, const t_DataType &THIRD, const t_DataType &FOURTH) :
        m_First( FIRST),
        m_Second( SECOND),
        m_Third( THIRD),
        m_Fourth( FOURTH)
      {
      }

      //! copy constructor
      VectorND( const VectorND &A) :
        m_First( A.m_First),
        m_Second( A.m_Second),
        m_Third( A.m_Third),
        m_Fourth( A.m_Fourth)
      {
      }

      //! move
      VectorND( VectorND && A) :
        m_First( std::move( A.m_First)),
        m_Second( std::move( A.m_Second)),
        m_Third( std::move( A.m_Third)),
        m_Fourth( std::move( A.m_Fourth))
      {
      }

      //! virtual copy constructor
      VectorND< 4, t_DataType> *Clone() const
      {
        return new VectorND< 4, t_DataType>( *this);
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

      //! return const refernce to first element
      t_DataType const &First() const
      {
        return m_First;
      }

      //! return reference to changable element
      t_DataType       &First()
      {
        return m_First;
      }

      //! return const refernce to second element
      t_DataType const &Second() const
      {
        return m_Second;
      }

      //! return reference to second element
      t_DataType       &Second()
      {
        return m_Second;
      }

      //! return const refernce to third element
      t_DataType const &Third() const
      {
        return m_Third;
      }

      //! return reference to third element
      t_DataType       &Third()
      {
        return m_Third;
      }

      //! return const refernce to fourth element
      t_DataType const &Fourth() const
      {
        return m_Fourth;
      }

      //! return reference to fourth element
      t_DataType       &Fourth()
      {
        return m_Fourth;
      }

    ///////////////
    // operators //
    ///////////////

      //! operator( POS) return const reference to element at POS
      t_DataType const &operator()( const size_t POS) const
      {
        switch( POS)
        {
          //first element
          case 0:
          {
            return m_First;
          }
          //second element
          case 1:
          {
            return m_Second;
          }
          //third element
          case 2:
          {
            return m_Third;
          }
          //fourth element
          case 3:
          {
            return m_Fourth;
          }
          //if pos is higher than number elements
          default:
          {
            BCL_Exit( "try to access non existing element in VectorND< 4, t_DataType>! " + util::Format()( POS) + " >= 4", -1);
            return m_First;
          }
        }
      }

      //! operator( POS) return reference to changeable element at POS
      t_DataType       &operator()( const size_t POS)
      {
        switch( POS)
        {
          //first element
          case 0:
          {
            return m_First;
          }
          //second element
          case 1:
          {
            return m_Second;
          }
          //third element
          case 2:
          {
            return m_Third;
          }
          //fourth element
          case 3:
          {
            return m_Fourth;
          }
          //if pos is higher than number elements
          default:
          {
            BCL_Exit( "try to access non existing element in VectorND< 4, t_DataType>! " + util::Format()( POS) + " >= 4", -1);
            return m_First;
          }
        }
      }

      //! @brief equal operator
      VectorND &operator =( const VectorND &VECTOR)
      {
        if( this != &VECTOR)
        {
          m_First = VECTOR.m_First;
          m_Second = VECTOR.m_Second;
          m_Third = VECTOR.m_Third;
          m_Fourth = VECTOR.m_Fourth;
        }
        return *this;
      }

      //! @brief move assignment operator
      VectorND &operator =( VectorND && VECTOR)
      {
        if( this != &VECTOR)
        {
          m_First = std::move( VECTOR.m_First);
          m_Second = std::move( VECTOR.m_Second);
          m_Third = std::move( VECTOR.m_Third);
          m_Fourth = std::move( VECTOR.m_Fourth);
        }
        return *this;
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! write VectorND< 4, t_DataType> to std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        //write data
        io::Serialize::Write( m_First, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Second, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Third, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Fourth, OSTREAM, INDENT);

        //return
        return OSTREAM;
      }

      //! read VectorND< 4, t_DataType> from std::istream
      std::istream &Read( std::istream &ISTREAM)
      {
        //read data
        io::Serialize::Read( m_First, ISTREAM);
        io::Serialize::Read( m_Second, ISTREAM);
        io::Serialize::Read( m_Third, ISTREAM);
        io::Serialize::Read( m_Fourth, ISTREAM);

        //end
        return ISTREAM;
      }

    }; // template class VectorND< 4, t_DataType>

    // instantiate s_Instance
    template< unsigned int t_Size, typename t_DataType>
    const util::SiPtr< const util::ObjectInterface> VectorND< t_Size, t_DataType>::s_Instance
    (
      GetObjectInstances().AddInstance( new VectorND< t_Size, t_DataType>())
    );

    // instantiate s_Instance
    template< typename t_DataType>
    const util::SiPtr< const util::ObjectInterface> VectorND< 1, t_DataType>::s_Instance
    (
      GetObjectInstances().AddInstance( new VectorND< 1, t_DataType>())
    );

    // instantiate s_Instance
    template< typename t_DataType>
    const util::SiPtr< const util::ObjectInterface> VectorND< 2, t_DataType>::s_Instance
    (
      GetObjectInstances().AddInstance( new VectorND< 2, t_DataType>())
    );

    // instantiate s_Instance
    template< typename t_DataType>
    const util::SiPtr< const util::ObjectInterface> VectorND< 3, t_DataType>::s_Instance
    (
      GetObjectInstances().AddInstance( new VectorND< 3, t_DataType>())
    );

    // instantiate s_Instance
    template< typename t_DataType>
    const util::SiPtr< const util::ObjectInterface> VectorND< 4, t_DataType>::s_Instance
    (
      GetObjectInstances().AddInstance( new VectorND< 4, t_DataType>())
    );

    //! this is a dummy *  operator to allow usage of SumFunction class with std::pair< AA, AA>
    template< typename T1, unsigned int N, typename t_DataType>
    inline
    VectorND< N, t_DataType>
    operator *( const T1 &VALUE, const VectorND< N, t_DataType> &VECTOR_ND)
    {
      BCL_Exit( "This dummy function should never be called!", -1);
      return VECTOR_ND;
    }

    //! this is a dummy += operator to allow usage of SumFunction class with std::pair< AA, AA>
    template< unsigned int N, typename t_DataType>
    inline
    VectorND< N, t_DataType>
    operator +=( const VectorND< N, t_DataType> &VECTOR_ND_A, const VectorND< N, t_DataType> &VECTOR_ND_B)
    {
      BCL_Exit( "This dummy function should never be called!", -1);
      return VECTOR_ND_A;
    }

    //! @brief operator == to compare to vector nd's of same size
    //! @param VECTOR_ND_LHS left hand side operand
    //! @param VECTOR_ND_RHS right hand side operand
    //! @return true if all elements are equal, false otherwise
    template< unsigned int N, typename t_DataType>
    inline
    bool
    operator ==( const VectorND< N, t_DataType> &VECTOR_ND_LHS, const VectorND< N, t_DataType> &VECTOR_ND_RHS)
    {
      // check all pairs for equality
      for( unsigned int i( 0); i < N; ++i)
      {
        if( VECTOR_ND_LHS( i) != VECTOR_ND_RHS( i))
        {
          return false;
        }
      }

      // end - all elements are equal
      return true;
    }

    //! @brief operator != to compare to vector nd's of same size
    //! @param VECTOR_ND_LHS left hand side operand
    //! @param VECTOR_ND_RHS right hand side operand
    //! @return true if all elements are equal, false otherwise
    template< unsigned int N, typename t_DataType>
    inline
    bool
    operator !=( const VectorND< N, t_DataType> &VECTOR_ND_LHS, const VectorND< N, t_DataType> &VECTOR_ND_RHS)
    {
      // check all pairs for equality
      for( unsigned int i( 0); i < N; ++i)
      {
        if( VECTOR_ND_LHS( i) != VECTOR_ND_RHS( i))
        {
          return true;
        }
      }

      // end - all elements are equal
      return false;
    }

    //! @brief operator < to compare to vector nd's of same size
    //! @param VECTOR_ND_LHS left hand side operand
    //! @param VECTOR_ND_RHS right hand side operand
    //! @return true if the first position in VECTOR_ND_LHS LHS that differs from VECTOR_ND_RHS is less
    template< unsigned int N, typename t_DataType>
    inline
    bool
    operator <( const VectorND< N, t_DataType> &VECTOR_ND_LHS, const VectorND< N, t_DataType> &VECTOR_ND_RHS)
    {
      // compare the first non-equal elements according to <
      for( unsigned int i( 0); i < N; ++i)
      {
        if( VECTOR_ND_LHS( i) != VECTOR_ND_RHS( i))
        {
          return VECTOR_ND_LHS( i) < VECTOR_ND_RHS( i);
        }
      }

      // end - all elements are equal
      return false;
    }

  } // namespace storage
} // namespace bcl

#endif //BCL_STORAGE_VECTOR_ND_H_

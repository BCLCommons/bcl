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

#ifndef BCL_MATH_TENSOR_H_
#define BCL_MATH_TENSOR_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_matrix_const_reference.h"
#include "linal/bcl_linal_matrix_operations.h"
#include "linal/bcl_linal_matrix_reference.h"
#include "linal/bcl_linal_vector_operations.h"
#include "linal/bcl_linal_vector_reference.h"
#include "util/bcl_util_object_instances.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Tensor
    //! @brief This class is a template class for Tensors
    //!
    //! @see Bjarne Stroustrup, "The C++ Programming Language", chapter 11.6, pp 282-283
    //! @see Bjarne Stroustrup, "The C++ Programming Language", chapter 22.4, pp 662-685
    //!
    //! @see @link example_math_tensor.cpp @endlink
    //! @author woetzen
    //! @date 16.11.2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class Tensor :
      public util::SerializableInterface
    {

    private:

    //////////
    // data //
    //////////

      linal::Vector< t_DataType> m_Data; //!< Vector to hold the underlying data
      size_t m_NumberLayers;             //!< number layers
      size_t m_NumberRows;               //!< number rows
      size_t m_NumberCols;               //!< number columns

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! construct Tensor from optional number of layers, optional number of rows, optional number of columns, and optional single element
      Tensor( const size_t LAYERS = 0, const size_t ROWS = 0, const size_t COLS = 0, const t_DataType &VALUE = t_DataType()) :
        m_Data( LAYERS * ROWS * COLS, VALUE),
        m_NumberLayers( LAYERS),
        m_NumberRows( ROWS),
        m_NumberCols( COLS)
      {
      }

      //! construct TensorBase from number of layers, number of rows, number of columns and pointer on data
      Tensor( const size_t LAYERS, const size_t ROWS, const size_t COLS, const t_DataType *DATA) :
        m_Data( LAYERS * ROWS * COLS, DATA),
        m_NumberLayers( LAYERS),
        m_NumberRows( ROWS),
        m_NumberCols( COLS)
      {
      }

      //! copy constructor
      Tensor( const Tensor< t_DataType> &PARENT) :
        m_Data( PARENT.m_Data),
        m_NumberLayers( PARENT.m_NumberLayers),
        m_NumberRows( PARENT.m_NumberRows),
        m_NumberCols( PARENT.m_NumberCols)
      {
      }

      //! virtual copy constructor
      virtual Tensor< t_DataType> *Clone() const
      {
        return new Tensor< t_DataType>( *this);
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

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const
      {
        static const std::string s_name( "tensor");
        return s_name;
      }

      //! return number of layers
      virtual size_t NumberLayers() const
      {
        return m_NumberLayers;
      }

      //! return number of rows
      virtual size_t GetNumberRows() const
      {
        return m_NumberRows;
      }

      //! return number of columns
      virtual size_t GetNumberCols() const
      {
        return m_NumberCols;
      }

      //! return pointer on begin
      virtual       t_DataType *Begin()
      {
        return m_Data.Begin();
      }

      //! return const pointer on begin
      virtual const t_DataType *Begin() const
      {
        return m_Data.Begin();
      }

      //! return pointer on end
      virtual       t_DataType *End()
      {
        return m_Data.End();
      }

      //! return const pointer on end
      virtual const t_DataType *End() const
      {
        return m_Data.End();
      }

      //! return the internal vector
      linal::VectorConstReference< t_DataType> GetValues() const
      {
        return linal::VectorConstReference< t_DataType>( m_Data);
      }

      //! return the internal vector
      linal::VectorReference< t_DataType> GetValues()
      {
        return linal::VectorReference< t_DataType>( m_Data);
      }

      //! C-style data access with [] gives a pointer on a LAYER
      t_DataType *operator[]( const size_t LAYER)
      {
        IsValidPosition( LAYER, 0, 0);
        return m_Data.Begin() + LAYER * m_NumberRows * m_NumberCols;
      }

      //! C-style data access with [] gives a const pointer on a LAYER
      const t_DataType *operator[]( const size_t LAYER) const
      {
        IsValidPosition( LAYER, 0, 0);
        return m_Data.Begin() + LAYER * m_NumberRows * m_NumberCols;
      }

      //! return reference to changeable element ( LAYER, ROW, COL)
            t_DataType &operator()( const size_t LAYER, const size_t ROW, const size_t COL)
      {
        return m_Data( ( LAYER * m_NumberRows + ROW) * m_NumberCols + COL);
      }

      //! return of element element ( LAYER, ROW, COL)
      const t_DataType &operator()( const size_t LAYER, const size_t ROW, const size_t COL) const
      {
        return m_Data( ( LAYER * m_NumberRows + ROW) * m_NumberCols + COL);
      }

      //! get a particular layer
      linal::MatrixConstReference< t_DataType> operator()( const size_t &LAYER) const
      {
        return linal::MatrixConstReference< t_DataType>( m_NumberRows, m_NumberCols, operator[]( LAYER));
      }

      //! get a particular layer
      linal::MatrixReference< t_DataType> operator()( const size_t &LAYER)
      {
        return linal::MatrixReference< t_DataType>( m_NumberRows, m_NumberCols, operator[]( LAYER));
      }

      //! return number of elements
      virtual size_t GetSize() const
      {
        return m_Data.GetSize();
      }

      //! return whether the tensor is empty
      virtual bool IsEmpty() const
      {
        return m_Data.IsEmpty();
      }

    //////////////////
    // manipulation //
    //////////////////

      //! set all elements in tensor to t_DataType( 0)
      virtual Tensor< t_DataType> &SetZero()
      {
        m_Data = t_DataType( 0);
        return *this;
      }

      //! set all elements in tensor to VALUE
      virtual Tensor< t_DataType> &SetValue( const t_DataType &VALUE)
      {
        m_Data = VALUE;
        return *this;
      }

      //! @brief normalize the tensor
      virtual Tensor< t_DataType> &Normalize()
      {
        m_Data.Normalize();

        // end
        return *this;
      }

    //////////////////
    // manipulation //
    //////////////////

      //! copies elements of argument MATRIX into this object at position ( LAYER)
      virtual Tensor< t_DataType> &ReplaceLayer( const size_t LAYER, const linal::MatrixConstInterface< t_DataType> &MATRIX);

    ////////////////
    // operations //
    ////////////////

      //! check whether tensor is qubic
      virtual bool IsQubic() const
      {
        return ( m_NumberLayers == m_NumberRows && m_NumberRows == m_NumberCols);
      }

      //! check whether tensor is square tensor of dimension NxN
      virtual bool IsTensorND( const size_t N) const
      {
        return ( IsQubic() && m_NumberLayers == N);
      }

      //! @brief create sub tensor from given position and width
      Tensor< t_DataType> SubTensor
      (
        const size_t NUM_LAYER, const size_t NUM_ROW, const size_t NUM_COL,
        const size_t LAYER = 0, const size_t ROW = 0, const size_t COL = 0
      ) const
      {
        // calculate last index
        const size_t last_layer( LAYER + NUM_LAYER);
        const size_t last_row(     ROW + NUM_ROW);
        const size_t last_col(     COL + NUM_COL);

        // check if the sub tensor is not too big
        IsValidPosition( LAYER, ROW, COL);
        IsValidPosition( last_layer, last_row, last_col);

        // create the sub tensor
        Tensor< t_DataType> sub_tensor( NUM_LAYER, NUM_ROW, NUM_COL, t_DataType( 0));
        t_DataType *ptr( sub_tensor.m_Data.Begin());

        // size of layers
        const size_t this_layer_size( m_NumberRows * m_NumberCols);

        // iterate over layer
        for
        (
          const t_DataType
            *dat_layer( m_Data.Begin() + LAYER * this_layer_size),
            *dat_layer_end( m_Data.Begin() + ( last_layer) * this_layer_size);
          dat_layer != dat_layer_end;
          dat_layer += this_layer_size
        )
        {
          // iterate over rows
          for
          (
            const t_DataType
              *dat_row( dat_layer + ROW * m_NumberCols),
              *dat_row_end( dat_layer + ( last_row) * m_NumberCols);
            dat_row != dat_row_end;
            dat_row += m_NumberCols
          )
          {
            std::copy( dat_row + COL, dat_row + last_col, ptr);
            ptr += NUM_COL;
          }
        }

        // end
        return sub_tensor;
      }

      //! @brief create a padded tensor
      //! @param LAYER_PADDING the number of additional layers to pad the tensor
      //! @param ROW_PADDING the number of additional rows to pad the tensor
      //! @param COL_PADDING the number of additional columns to pad the tensor
      //! @param FILL_VALUE addign every pad element that value
      //! @return Tensor< t_DataType> that has additional layers, rows and cols preserving the tensor in the upper left corner
      Tensor< t_DataType> CreatePaddedTensor
      (
        const size_t LAYER_PADDING,
        const size_t ROW_PADDING,
        const size_t COL_PADDING,
        const t_DataType &FILL_VALUE = t_DataType( 0)
      ) const
      {
        const size_t pad_num_rows( m_NumberRows + ROW_PADDING);
        const size_t pad_num_cols( m_NumberCols + COL_PADDING);

        // create tensor
        Tensor< t_DataType> padded_tensor( m_NumberLayers + LAYER_PADDING, pad_num_rows, pad_num_cols, FILL_VALUE);
        t_DataType *ltarget( padded_tensor.m_Data.Begin());

        // copy layer by layer
        for
        (
          const t_DataType *ldata( m_Data.Begin()), *ldata_end( m_Data.End());
          ldata != ldata_end;
          ldata += m_NumberRows * m_NumberCols, ltarget += pad_num_rows * pad_num_cols
        )
        {
          t_DataType *rtarget( ltarget);
          // copy row by row
          for
          (
            const t_DataType *rdata( ldata), *rdata_end( ldata + m_NumberRows * m_NumberCols);
            rdata != rdata_end;
            rdata += m_NumberCols, rtarget += pad_num_cols
          )
          {
            std::copy
            (
              rdata,
              rdata + m_NumberCols,
              rtarget
            );
          }
        }

        // end
        return padded_tensor;
      }

    ///////////////
    // operators //
    ///////////////

      //! set all elements of tensor to one given VALUE
      virtual Tensor< t_DataType> &operator =( const t_DataType &VALUE)
      {
        m_Data =( VALUE);
        return *this;
      }

      //! operator += VALUE
      virtual Tensor< t_DataType> &operator +=( const t_DataType &VALUE)
      {
        m_Data += ( VALUE);
        return *this;
      }

      //! operator -= VALUE
      virtual Tensor< t_DataType> &operator -=( const t_DataType &VALUE)
      {
        m_Data -= ( VALUE);
        return *this;
      }

      //! operator *= VALUE
      virtual Tensor< t_DataType> &operator *=( const t_DataType &VALUE)
      {
        m_Data *= ( VALUE);
        return *this;
      }

      //! operator /= VALUE
      virtual Tensor< t_DataType> &operator /=( const t_DataType &VALUE)
      {
        m_Data /= ( VALUE);
        return *this;
      }

      //! copy one tensor into another tensor
      virtual Tensor< t_DataType> &operator =( const Tensor< t_DataType> &TENSOR)
      {
        //copy the dimensions fo the tensor
        m_NumberLayers = TENSOR.NumberLayers();
        m_NumberRows   = TENSOR.GetNumberRows();
        m_NumberCols   = TENSOR.GetNumberCols();

        //call base class operator
        m_Data = TENSOR.m_Data;
        return *this;
      }

      //! @brief logical equal operator
      //! @param TENSOR_RHS right hand side tensor
      //! @return true if the dimensions and all values match
      bool operator ==( const Tensor< t_DataType> &TENSOR_RHS) const
      {
        return
        (
             m_NumberCols   == TENSOR_RHS.m_NumberCols
          && m_NumberRows   == TENSOR_RHS.m_NumberRows
          && m_NumberLayers == TENSOR_RHS.m_NumberLayers
          && m_Data         == TENSOR_RHS.m_Data
        );
      }

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
        {
          io::Serializer serializer;
          serializer.SetClassDescription( "template class for tensors");
          serializer.AddInitializer
            (
             "data",
             "vector to hold the underlying data",
             io::Serialization::GetAgent( &m_Data)
             );
          serializer.AddInitializer
            (
             "number layers",
             "number of layers",
             io::Serialization::GetAgent( &m_NumberLayers)
             );
          serializer.AddInitializer
            (
             "number rows",
             "number of rows",
             io::Serialization::GetAgent( &m_NumberRows)
             );
          serializer.AddInitializer
            (
             "number columns",
             "number of columns",
             io::Serialization::GetAgent( &m_NumberCols)
             );

          return serializer;
      }
    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! write tensor to std::ostream using the given util::Format
      virtual std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! read tensor from io::IFStream
      virtual std::istream &Read( std::istream &ISTREAM);

      //! check whether position is valid
      void IsValidPosition( const size_t LAYER, const size_t ROW, const size_t COL) const
      {
        BCL_Assert( LAYER < m_NumberLayers, "LAYER extend size of Matrix " + util::Format()( LAYER) + " > [0.." + util::Format()( m_NumberLayers) + "]!");
        BCL_Assert( ROW   < m_NumberRows  ,   "ROW extend size of Matrix " + util::Format()(   ROW) + " > [0.." + util::Format()( m_NumberRows  ) + "]!");
        BCL_Assert( COL   < m_NumberCols  ,   "COL extend size of Matrix " + util::Format()(   COL) + " > [0.." + util::Format()( m_NumberCols  ) + "]!");
      }

    }; //end template Tensor class

    // instantiate s_Instance
    template< typename t_DataType>
    const util::SiPtr< const util::ObjectInterface> Tensor< t_DataType>::s_Instance
    (
      GetObjectInstances().AddInstance( new Tensor< t_DataType>())
    );

    //! copies elements of argument MATRIX into this object at position (LAYER)
    template< typename t_DataType> Tensor< t_DataType> &Tensor< t_DataType>::ReplaceLayer( const size_t LAYER, const linal::MatrixConstInterface< t_DataType> &MATRIX)
    {
      IsValidPosition( LAYER, 0, 0);
      BCL_Assert( MATRIX.GetNumberRows() == m_NumberRows && MATRIX.GetNumberCols() == m_NumberCols, "given Layer has wrong dimensions!");

      std::copy
      (
        MATRIX.Begin(),
        MATRIX.End(),
        m_Data.Begin() + LAYER * m_NumberRows * m_NumberCols
      );

      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! write Tensor to STREAM using the given FORMAT
    template< typename t_DataType> inline std::ostream &Tensor< t_DataType>::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write header
      io::Serialize::Write( m_NumberLayers, OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_NumberRows, OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_NumberCols, OSTREAM, INDENT) << '\n';

      // write elements
      for( size_t i( 0); i < NumberLayers(); ++i)
      {
        for( size_t j( 0); j < GetNumberRows(); ++j)
        {
          for( size_t k( 0); k < GetNumberCols(); ++k)
          {
            io::Serialize::Write( operator()( i, j, k), OSTREAM, INDENT) << '\t';
          }
          OSTREAM << '\n';
        }
        OSTREAM << '\n';
      }

      // end
      return OSTREAM;
    }

    //! read Tensor from std::istream
    template< typename t_DataType> inline std::istream &Tensor< t_DataType>::Read( std::istream &ISTREAM)
    {
      // read data
      io::Serialize::Read( m_NumberLayers, ISTREAM);
      io::Serialize::Read( m_NumberRows, ISTREAM);
      io::Serialize::Read( m_NumberCols, ISTREAM);
      m_Data = linal::Vector< t_DataType>( m_NumberLayers * m_NumberRows * m_NumberCols);

      for( size_t i( 0); i < m_NumberLayers; ++i)
      {
        for( size_t j( 0); j < m_NumberRows; ++j)
        {
          for( size_t k( 0); k < m_NumberCols; ++k)
          {
            io::Serialize::Read( operator()( i, j, k), ISTREAM);
          }
        }
      }

      // return
      return ISTREAM;
    }

    BCL_EXPIMP_TEMPLATE template class BCL_API Tensor< double>;

  } // namespace math
} // namespace bcl

#endif //BCL_MATH_TENSOR_H_

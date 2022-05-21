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

#ifndef BCL_UTIL_VOXEL_GRID_H_
#define BCL_UTIL_VOXEL_GRID_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_logger_interface.h"
#include "bcl_util_si_ptr_vector.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_running_min_max.h"
#include "math/bcl_math_statistics.h"
#include "storage/bcl_storage_triplet.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class VoxelGrid
    //! @brief Optimal space hashing class for rapid identification of neighbors
    //! @details Provides retrieval of objects in 3D space by space hashing (1,2, or 3D Voxel grid depending on object
    //!          density. Object references and coordinates are stored on the grid and allow for very rapid retrieval of
    //!          neighbors within a distance <= the resolution of the grid.
    //! @see @link example_score_slicelist_benchmark.cpp @endlink
    //! @author mendenjl
    //! @date Nov 19, 2016
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_DataType>
    class VoxelGrid :
      public ObjectInterface
    {

    private:

      enum BorderFlag
      {
        e_None = 0,
        e_XMin = 1,
        e_XMax = 2,
        e_YMin = 4,
        e_YMax = 8,
        e_ZMin = 16,
        e_ZMax = 32
      };

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class ObjReference
      //! @brief Reference to each object and its associated coordinates. We could use a pair instead, but item names
      //!        are more descriptive
      //! @author mendenjl
      //! @date Nov 29, 2016
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////
      class ObjReference :
        public ObjectInterface
      {
      public:

      //////////
      // data //
      //////////

        SiPtr< const t_DataType>      m_Object;
        SiPtr< const linal::Vector3D> m_Coordinates;

      public:

        //! @brief default constructor
        ObjReference
        (
          const SiPtr< const t_DataType> &OBJ = SiPtr< const t_DataType>(),
          const SiPtr< const linal::Vector3D> &COORD = SiPtr< const linal::Vector3D>()
        ) :
          m_Object( OBJ),
          m_Coordinates( COORD)
        {
        }

        //! virtual copy constructor
        ObjReference *Clone() const
        {
          return new ObjReference( *this);
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

        //! @brief true if the object is defined
        bool IsDefined() const
        {
          return m_Coordinates.IsDefined();
        }

        //! read from std::istream ( dummy)
        std::istream &Read( std::istream &ISTREAM)
        {
          BCL_ExitWithoutCallstack( "Not implemented", -1);
          return ISTREAM;
        }

        //! write to std::ostream
        std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
        {
          BCL_ExitWithoutCallstack( "Not implemented", -1);
          return OSTREAM;
        }

      };

    //////////
    // data //
    //////////

      storage::Vector< storage::Vector< ObjReference> > m_Assignments; //!< Object Assignments

      //! neighboring voxels; used only for 2-3D grids, ignored for 1D. This set is for directed edges
      //! used when detected internal clashes
      storage::Vector< SiPtrVector< const storage::Vector< ObjReference> > > m_Edges;

      //! neighboring voxels; used only for 2-3D grids, ignored for 1D
      //! This set is for undirected edges used when detected external clashes
      storage::Vector< SiPtrVector< const storage::Vector< ObjReference> > > m_EdgesComplete;

      //! Border flags, used if not caching edges
      std::vector< int> m_BorderFlags;

      size_t m_NBinsX;        //!< Number of bins each direction
      size_t m_NBinsY;
      size_t m_NBinsZ;
      double m_Resolution;    //!< Requested resolution by the user
      bool   m_CacheEdges;    //!< Should be true if the same grid will be used repeatedly (~20% speed improvement)
                              //!  Edges should not be cached if the grid is only going to be used once, otherwise a
                              //!  significant speed penalty can be expected
      double m_ResolutionX;   //!< Resolution of the optimized grid in the X direction
      double m_ResolutionY;   //!< Resolution of the optimized grid in the Y direction
      double m_ResolutionZ;   //!< Resolution of the optimized grid in the Z direction
      double m_MinResolution; //!< Minimum of resolution in any axes
      double m_MinX;          //!< Grid origin
      double m_MinY;
      double m_MinZ;
      math::RunningMinMax< linal::Vector3D> m_Box; //!< Bounding box for points given
      size_t m_NDimensional; //!< Number of dimensions expanded by the voxel grid
      size_t m_NItems;       //!< Number of items
      size_t m_MinNumberElements; //!< Minimum number of elements this object should have; if smaller, reallocate grid
      size_t m_MaxNumberElements; //!< Maximum number of elements this object should have; if smaller, reallocate grid

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      //! @param RESOLUTION resolution of the grid
      //! @param CACHE_EDGES true to cache edges between gridpoints for 2/3D voxel grids
      VoxelGrid( const double &RESOLUTION = 4.0, const bool &CACHE_EDGES = false) :
        m_NBinsX( 0),
        m_NBinsY( 0),
        m_NBinsZ( 0),
        m_Resolution( RESOLUTION),
        m_CacheEdges( CACHE_EDGES),
        m_ResolutionX( RESOLUTION),
        m_ResolutionY( RESOLUTION),
        m_ResolutionZ( RESOLUTION),
        m_MinResolution( RESOLUTION),
        m_MinX( 0),
        m_MinY( 0),
        m_MinZ( 0),
        m_NDimensional( 0),
        m_NItems( 0),
        m_MinNumberElements( 0),
        m_MaxNumberElements( 0)
      {
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief extract the 3D coordinate of a given t_DataType input. TO BE IMPLEMENTED BY DERIVED CLASSES.
      //! @param INPUT pointer to the t_DataType
      //! @return 3D Vector of the input's coordinates
      virtual SiPtr< const linal::Vector3D> ExtractPosition( const SiPtr< const t_DataType> &INPUT) const = 0;

      //! @brief extract the 3D coordinates of a given t_DataType input. TO BE IMPLEMENTED BY DERIVED CLASSES if one input can have multiple positions
      //! @param INPUT pointer to the t_DataType
      //! @return 3D Vector of the input's coordinates
      virtual SiPtrVector< const linal::Vector3D> ExtractPositions( const SiPtr< const t_DataType> &INPUT) const
      {
        auto siptr( ExtractPosition( INPUT));
        return
          siptr.IsDefined() ? SiPtrVector< const linal::Vector3D>( size_t( 1), &*siptr) : SiPtrVector< const linal::Vector3D>();
      }

      //! @brief check if two list items are the same. TO BE IMPLEMENTED BY DERIVED CLASSES.
      //! @note for most cases, this can just return the == comparison value, but for some objects,
      //!       this might not be sufficient (e.g. AAs have to be compared via ChainIds and SeqIDs)
      //! @param ITEM_1 first comparison item
      //! @param ITEM_2 second comparison item
      //! @return true if list items are the same, false otherwise
      virtual bool IsSameItem( const t_DataType &ITEM_1, const t_DataType &ITEM_2) const = 0;

      //! @brief Get the number of items in the voxel grid
      size_t GetNumberItems() const
      {
        return m_NItems;
      }

      //! @brief Get the actual internally used dimension of the voxel grid
      size_t GetDimension() const
      {
        return m_NDimensional;
      }

      //! @brief Find neighbors within a specified distance of a given DATA
      //! @param INPUT datapoint of interest
      //! @param NEIGHBORHOOD, must be <= resolution
      storage::Vector< storage::Pair< SiPtr< const t_DataType>, double> > GetNeighbors( const t_DataType &INPUT, const double &NEIGHBORHOOD)
      {
        return m_NDimensional <= size_t( 1) ? GetNeighbors1D( INPUT, NEIGHBORHOOD) : GetNeighborsMultiDimensional( INPUT, NEIGHBORHOOD);
      }

      //! @brief Find neighbors within a specified distance of a given DATA
      //! @param INPUT datapoint of interest
      //! @param NEIGHBORHOOD, must be <= resolution
      storage::Vector< storage::Pair< SiPtr< const t_DataType>, double> > GetNeighbors( const linal::Vector3D &INPUT, const double &NEIGHBORHOOD)
      {
        return m_NDimensional <= size_t( 1) ? GetNeighbors1D( INPUT, NEIGHBORHOOD) : GetNeighborsMultiDimensional( INPUT, NEIGHBORHOOD);
      }

      //! @brief Find all neighbor pairs within a specified distance
      //! @param GRID other grid to test for neighbors between
      //! @param NEIGHBORHOOD, must be <= resolution
      storage::Vector< storage::Triplet< SiPtr< const t_DataType>, SiPtr< const t_DataType>, double> >
        GetNeighborsIn( const VoxelGrid &GRID, const double &NEIGHBORHOOD) const
      {
        return m_NDimensional <= size_t( 1) ? GetNeighbors1D( GRID, NEIGHBORHOOD) : GetNeighborsMultiDimensional( GRID, NEIGHBORHOOD);
      }

      //! @brief Find all neighbor pairs within a specified distance
      //! @param NEIGHBORHOOD, must be <= resolution
      storage::Vector< storage::Triplet< SiPtr< const t_DataType>, SiPtr< const t_DataType>, double> >
        GetNeighbors( const double &NEIGHBORHOOD) const
      {
        return m_NDimensional <= size_t( 1) ? GetNeighbors1D( NEIGHBORHOOD) : GetNeighborsMultiDimensional( NEIGHBORHOOD);
      }

      //! @brief Update the VoxelGrid with new data
      //! @param NEW_DATA SiPtrVector of the new data
      virtual void SetObjects( const SiPtrVector< const t_DataType> &NEW_DATA)
      {
        // track the boundaries of the cube enclosing all points given
        m_Box.Reset();

        // store references to the objects and their coordinates: on AAs getting the position of the AA
        // is slow so not caching this here results in a substantial slowdown (>20% often)
        storage::Vector< ObjReference> ref_vec;
        ref_vec.AllocateMemory( NEW_DATA.GetSize());
        for( auto itr( NEW_DATA.Begin()), itr_end( NEW_DATA.End()); itr != itr_end; ++itr)
        {
          SiPtrVector< const linal::Vector3D> coord( ExtractPositions( **itr));
          for( auto itr_c( coord.Begin()), itr_c_end( coord.End()); itr_c != itr_c_end; ++itr_c)
          {
            m_Box += **itr_c;
            ref_vec.PushBack( ObjReference( *itr, **itr_c));
          }
        }

        // Determine optimal dimensionality for the grid and set it up
        SetupGrid( m_Box.GetMin(), m_Box.GetMax(), ref_vec.GetSize());

        // populate the grid
        for( auto itr( ref_vec.Begin()), itr_end( ref_vec.End()); itr != itr_end; ++itr)
        {
          const size_t index( GetIndex( *itr->m_Coordinates));
          m_Assignments( index).PushBack( *itr);
        }
        m_NItems = ref_vec.GetSize();
      }

    ///////////////////////
    // data manipulation //
    ///////////////////////

      //! @brief insert a new t_DataType object into the manager
      //! @param NEW_ITEM pointer to item to insert
      void InsertObject( const t_DataType &NEW_ITEM)
      {
        SiPtr< const linal::Vector3D> coord( ExtractPosition( NEW_ITEM));
        if( coord.IsDefined())
        {
          const size_t index( GetIndex( *coord));
          m_Assignments( index).PushBack( ObjReference( NEW_ITEM, coord));
          ++m_NItems;
          m_Box += *coord;
        }
      }

      //! @brief inserts multiple t_DataType objects into the manager
      //! @param NEW_ITEMS pointer vector to items to insert
      void InsertObjects( const SiPtrVector< const t_DataType> NEW_ITEMS)
      {
        for( size_t i( 0), sz( NEW_ITEMS.GetSize()); i < sz; ++i)
        {
          InsertObject( *NEW_ITEMS( i));
        }
      }

      //! @brief remove a t_DataType object from the manager
      //! @param ITEM_TO_REMOVE pointer to the item to remove
      void RemoveObject( const t_DataType &ITEM_TO_REMOVE)
      {
        SiPtr< const linal::Vector3D> pos( ExtractPosition( ITEM_TO_REMOVE));
        if( pos.IsDefined() && m_NBinsX)
        {
          const size_t index( GetIndex( *pos));
          storage::Vector< ObjReference> &vec( m_Assignments( index));

          for( size_t i( 0), sz( vec.GetSize()); i < sz; ++i)
          {
            if( IsSameItem( *vec( i).m_Object, ITEM_TO_REMOVE))
            {
              vec.RemoveElements( i, 1);
              --m_NItems;
              break;
            }
          }
        }
      }

      //! @brief remove multiple t_DataType objects from the manager
      //! @param ITEMS_TO_REMOVE pointer vector to items to remove
      void RemoveObjects( const SiPtrVector< const t_DataType> &ITEMS_TO_REMOVE)
      {
        for( size_t i( 0), sz( ITEMS_TO_REMOVE.GetSize()); i < sz; ++i)
        {
          RemoveObject( *ITEMS_TO_REMOVE( i));
        }
      }

      //! @brief Remove all objects from the voxel grid while keeping the grid itself intact
      void Clear()
      {
        for( size_t i( 0), sz( m_Assignments.GetSize()); i < sz; ++i)
        {
          m_Assignments( i).Reset();
        }
        m_NItems = 0;
      }

      //! @brief Translate the grid
      //! @param TRANSLATION amount to translate by
      void Translate( const linal::Vector3D &TRANSLATION)
      {
        if( m_NItems)
        {
          m_MinX += TRANSLATION.X();
          m_MinY += TRANSLATION.Y();
          m_MinZ += TRANSLATION.Z();
          linal::Vector3D new_min( m_Box.GetMin());
          new_min += TRANSLATION;
          linal::Vector3D new_max( m_Box.GetMax());
          new_max += TRANSLATION;
          m_Box.Reset();
          m_Box += new_min;
          m_Box += new_max;
        }
      }

    ////////////////////
    // input & output //
    ////////////////////

      //! @brief GetClassIdentifier returns class name of the object
      //! @return returns string with the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    private:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

      //! @brief Get the position for the given coordinate vector. Places out-of-bounds vector in the bin closest to
      //!        where they belong
      //! @param X, Y, Z the position of interest
      //! @return index in m_Coordinates vector for the given point
      size_t GetIndex( const double &X, const double &Y, const double &Z) const
      {
        return size_t( std::max( std::min( ( X - m_MinX) / m_ResolutionX, m_NBinsX - 1.0), 0.0)) // x bin
               + m_NBinsX *
               (
                 size_t( std::max( std::min( ( Y - m_MinY) / m_ResolutionY, m_NBinsY - 1.0), 0.0))   // y bin
                 + size_t( std::max( std::min( ( Z - m_MinZ) / m_ResolutionZ, m_NBinsZ - 1.0), 0.0)) * m_NBinsY // z bin
               );
      }

      //! @brief Get the position for the given coordinate vector. Places out-of-bounds vector in the bin closest to
      //!        where they belong
      //! @param POS the position of interest
      //! @return index in m_Coordinates vector for the given point
      size_t GetIndex( const linal::Vector3D &POS) const
      {
        return this->GetIndex( POS.X(), POS.Y(), POS.Z());
      }

      //! @brief Try to allocate a grid of a given size
      //! If this size would cause the object to definitely exceed the max memory allotment, the grid will be truncated
      //! to include only the core. It is assumed that m_Resolution has been setup prior to this
      void SetupGrid( const linal::Vector3D &MIN_VALS, const linal::Vector3D &MAX_VALS, const size_t &EXPECTED_N_ELEMENTS)
      {
        // always clear the grid of any objects
        Clear();

        if( !EXPECTED_N_ELEMENTS)
        {
          return;
        }

        BCL_Assert
        (
          MIN_VALS( 0) <= MAX_VALS( 0) && MIN_VALS( 1) <= MAX_VALS( 1) && MIN_VALS( 2) <= MAX_VALS( 2),
          "Min vals and max vals incorrect"
        );

        // check whether we can continue using the existing grid
        if
        (
          MIN_VALS.X() >= m_MinX && MIN_VALS.Y() >= m_MinY && MIN_VALS.Z() >= m_MinZ
          && MAX_VALS.X() <= m_MinX + m_NBinsX * m_ResolutionX
          && MAX_VALS.Y() <= m_MinY + m_NBinsY * m_ResolutionY
          && MAX_VALS.Z() <= m_MinZ + m_NBinsZ * m_ResolutionZ
          && EXPECTED_N_ELEMENTS >= m_MinNumberElements
          && EXPECTED_N_ELEMENTS <= m_MaxNumberElements
        )
        {
          return;
        }

        linal::Vector3D diff( MAX_VALS - MIN_VALS);

        // Overallocation factor (>= 1.0). Grid dimensions overallocated by this amount to reduce allocation frequency
        static const double s_OverAlloc( m_CacheEdges ? 1.1 : 1.0);
        m_ResolutionX = m_ResolutionY = m_ResolutionZ = m_Resolution;

        m_NBinsX = std::max( size_t( ( diff.X() + m_ResolutionX) * s_OverAlloc / m_ResolutionX), size_t( 1));
        m_NBinsY = std::max( size_t( ( diff.Y() + m_ResolutionY) * s_OverAlloc / m_ResolutionY), size_t( 1));
        m_NBinsZ = std::max( size_t( ( diff.Z() + m_ResolutionZ) * s_OverAlloc / m_ResolutionZ), size_t( 1));

        // two bins are always less efficient than 1 - results in the same number of comparisons but a lot more iteration
        // In numerous tests, 3 bins is also virtually always much slower owing to the additional iterations and fact that
        // the central bin will still be compared with bins on both sides, so only minimal comparisons are avoided
        if( m_NBinsX <= 3)
        {
          m_NBinsX = 1;
          m_ResolutionX = std::numeric_limits< double>::max() / 3.0;
        }
        if( m_NBinsY <= 3)
        {
          m_NBinsY = 1;
          m_ResolutionY = std::numeric_limits< double>::max() / 3.0;
        }
        if( m_NBinsZ <= 3)
        {
          m_NBinsZ = 1;
          m_ResolutionZ = std::numeric_limits< double>::max() / 3.0;
        }

        // Next, decide whether this should be 1D, 2D, or 3D Voxel grid
        // Generally, the goal is to keep the expected objects per bin and bins per object as close to 1 as possible
        // There is also some bonus to keeping the dimensions smaller since it means vastly fewer bin iterations
        // Empirically, the benefit is roughly 15 to go from 2D - 1D (owing to a space-oversampling optimization that is
        // fast at 1D), and another 50 to go from 3D -> 2D. In other words, 3D clustering only makes sense when we have, on average,
        // more than about 65 elements total in each 1D grid of equivalent dimension.
        const double elements( EXPECTED_N_ELEMENTS);
        storage::Vector< double> dimension_bins( 7);
        storage::Vector< double> dimension_score( 7);
        double dimension_penalty[ 4] = { 0, 0, 15, 65};

        double best_score( elements + 2);
        int best_pos( 0);
        for( int pos( 1); pos < 8; ++pos)
        {
          const int use_x( !!( pos & 1)), use_y( !!( pos & 2)), use_z( !!( pos & 4));
          const int dimension( use_x + use_y + use_z);
          const double bins( ( use_x ? m_NBinsX : 1) * ( use_y ? m_NBinsY : 1) * ( use_z ? m_NBinsZ : 1));
          dimension_bins( pos - 1) = bins;
          // score is elements per bin + bins per element + dimension penalty
          const double score( elements / bins + bins / elements + dimension_penalty[ dimension]);
          dimension_score( pos - 1) = score;
          if( score < best_score)
          {
            best_pos = pos;
            best_score = score;
          }
        }

        if( !( best_pos & 1))
        {
          m_NBinsX = 1;
          m_ResolutionX = std::numeric_limits< double>::max() / 3.0;
        }
        if( !( best_pos & 2))
        {
          m_NBinsY = 1;
          m_ResolutionY = std::numeric_limits< double>::max() / 3.0;
        }
        if( !( best_pos & 4))
        {
          m_NBinsZ = 1;
          m_ResolutionZ = std::numeric_limits< double>::max() / 3.0;
        }
        m_NDimensional = ( ( m_NBinsX > 1) + ( m_NBinsY > 1) + ( m_NBinsZ > 1));

        // while it can be computed exactly when the score would be better for a different grid allocation, an
        // approximate method is adequate and better suited to ensuring that we avoid re-allocating the grid due to
        // changes in density as much as possible
        m_MinNumberElements = EXPECTED_N_ELEMENTS / 4;
        m_MaxNumberElements = m_NDimensional < 3 ? EXPECTED_N_ELEMENTS * 4 : std::numeric_limits< size_t>::max();

        // oversample if grid will be in 1D and there are sufficient elements.
        // This reduces the number of distance calculations and comparisons necessary
        if( m_NDimensional == size_t( 1))
        {
          // for 1D grids, it's not much of a burden to allocate for much larger grid than is absolutely required
          const double additional_over_alloc_factor( m_CacheEdges ? 1.25 : 1.0);
          if( m_NBinsX > size_t( 1))
          {
            m_NBinsX = size_t( double( m_NBinsX) * additional_over_alloc_factor);
            if( m_NBinsX * 3 < EXPECTED_N_ELEMENTS)
            {
              m_NBinsX *= 3;
              m_ResolutionX /= 3;
            }
          }
          else if( m_NBinsY > size_t( 1))
          {
            //  m_NBinsY = size_t( double( m_NBinsY) * additional_over_alloc_factor);
            if( m_NBinsY * 3 < EXPECTED_N_ELEMENTS)
            {
              m_NBinsY *= 3;
              m_ResolutionY /= 3;
            }
          }
          else // if( m_NBinsZ > size_t( 1))
          {
            //  m_NBinsZ = size_t( double( m_NBinsZ) * additional_over_alloc_factor);
            if( m_NBinsZ * 3 < EXPECTED_N_ELEMENTS)
            {
              m_NBinsZ *= 3;
              m_ResolutionZ /= 3;
            }
          }
        }
//        BCL_Debug( m_NBinsX);
//        BCL_Debug( m_NBinsY);
//        BCL_Debug( m_NBinsZ);
//        BCL_Debug( m_Box.GetMin());
//        BCL_Debug( m_Box.GetMax());
//        BCL_Debug( EXPECTED_N_ELEMENTS);

        const size_t total_requested_bins( m_NBinsX * m_NBinsY * m_NBinsZ);
        m_MinResolution = std::min( m_ResolutionX, std::min( m_ResolutionY, m_ResolutionZ));

        // define origin
        m_MinX = MIN_VALS( 0);
        m_MinY = MIN_VALS( 1);
        m_MinZ = MIN_VALS( 2);

        // update origin shifted (to more negative values) due to overallocation and dimensional collapsing
        m_MinX -= ( m_MinX + m_ResolutionX * m_NBinsX - MAX_VALS.X()) * 0.5;
        m_MinY -= ( m_MinY + m_ResolutionY * m_NBinsY - MAX_VALS.Y()) * 0.5;
        m_MinZ -= ( m_MinZ + m_ResolutionZ * m_NBinsZ - MAX_VALS.Z()) * 0.5;
        m_Assignments.Resize( total_requested_bins);

        if( m_NDimensional > size_t( 1))
        {
          // # of bins to explore in each direction. This can be easily expanded to try different binning strategies
          const int bins_to_x( 1), bins_to_y( 1), bins_to_z( 1);

          // last xyz bins, cached for convenience
          const int last_x_bin( m_NBinsX - 1), last_y_bin( m_NBinsY - 1), last_z_bin( m_NBinsZ - 1);

          // compute edges for 2-3 dimensional grids
          if( m_CacheEdges)
          {
            m_Edges.Resize( total_requested_bins);
            m_EdgesComplete.Resize( total_requested_bins);

            // distance in memory between adjacent bins on Z-axis
            const size_t z_stride( m_NBinsX * m_NBinsY);

            // determine side of edge vectors
            const size_t edge_vec_size_internal( m_NDimensional == size_t( 2) ? 4 : 13);
            const size_t edge_vec_size_external( m_NDimensional == size_t( 2) ? 8 : 26);

            // compute edges for multi-dimensional grids
            for( int z( 0), pos( 0), nbz( m_NBinsZ); z < nbz; ++z)
            {
              const int min_z( -std::min( z, bins_to_z)), max_z( std::min( z + bins_to_z, last_z_bin) - z + 1);
              for( int y( 0), nby( m_NBinsY); y < nby; ++y)
              {
                const int min_y( -std::min( y, bins_to_y)), max_y( std::min( y + bins_to_y, last_y_bin) - y + 1);
                for( int x( 0), nbx( m_NBinsX); x < nbx; ++x, ++pos)
                {
                  // determine extents in each direction
                  const int min_x( -std::min( x, bins_to_x)), max_x( std::min( x + bins_to_x, last_x_bin) - x + 1);

                  // reset existing edges
                  SiPtrVector< const storage::Vector< ObjReference> > &edges( m_Edges( pos));
                  edges.Reset();
                  edges.AllocateMemory( edge_vec_size_internal);
                  SiPtrVector< const storage::Vector< ObjReference> > &edgesc( m_EdgesComplete( pos));
                  edgesc.Reset();
                  edgesc.AllocateMemory( edge_vec_size_external);
                  for( int z_off( min_z), fz_off( min_z * z_stride); z_off < max_z; ++z_off, fz_off += z_stride)
                  {
                    for( int y_off( min_y), zy_off( fz_off + min_y * m_NBinsX); y_off < max_y; ++y_off, zy_off += m_NBinsX)
                    {
                      for( int x_off( min_x), zyx_off( zy_off + min_x); x_off < max_x; ++x_off, ++zyx_off)
                      {
                        if( zyx_off > 0)
                        {
                          // directed edges
                          edges.PushBack( ToSiPtr( m_Assignments( pos + zyx_off)));
                        }
                        if( zyx_off != 0)
                        {
                          // undirected edges; just ensure that it isn't the same point
                          edgesc.PushBack( ToSiPtr( m_Assignments( pos + zyx_off)));
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          else
          {
            m_BorderFlags.resize( total_requested_bins);
            // compute flags for multidimensional grids
            for( int z( 0), pos( 0), nbz( m_NBinsZ); z < nbz; ++z)
            {
              int z_flag( int( z ? e_None : e_ZMin) | int( z < last_z_bin ? e_None : e_ZMax));
              for( int y( 0), nby( m_NBinsY); y < nby; ++y)
              {
                int zy_flag( z_flag | int( y ? e_None : e_YMin) | int( y < last_y_bin ? e_None : e_YMax));
                for( int x( 0), nbx( m_NBinsX); x < nbx; ++x, ++pos)
                {
                  // determine extents in each direction
                  m_BorderFlags[ pos] = zy_flag | int( x ? e_None : e_XMin) | int( x < last_x_bin ? e_None : e_XMax);
                }
              }
            }
          }
        }
      }

      //! @brief Compute distance from this point to the nearest edge of the bounding box of all given points,
      //!        provided the point is outside the box
      //! @param X the coordinate of interest
      double GetSqDistanceOutsideBoundingBox( const linal::Vector3D &X) const
      {
        const linal::Vector3D &mn( m_Box.GetMin()), &mx( m_Box.GetMax());
        double oob_dist_sq( 0.0);
        if( X.X() < mn.X())
        {
          oob_dist_sq += math::Sqr( X.X() - mn.X());
        }
        else if( X.X() > mx.X())
        {
          oob_dist_sq += math::Sqr( X.X() - mx.X());
        }
        if( X.Y() < mn.Y())
        {
          oob_dist_sq += math::Sqr( X.Y() - mn.Y());
        }
        else if( X.Y() > mx.Y())
        {
          oob_dist_sq += math::Sqr( X.Y() - mx.Y());
        }
        if( X.Z() < mn.Z())
        {
          oob_dist_sq += math::Sqr( X.Z() - mn.Z());
        }
        else if( X.Z() > mx.Z())
        {
          oob_dist_sq += math::Sqr( X.Z() - mx.Z());
        }
        return oob_dist_sq;
      }

      //! @brief Compute distance from this point to the nearest edge of the bounding box of all given points,
      //!        provided the point is outside the box
      //! @param X the coordinate of interest
      double GetSqDistanceOutsideBoundingBox( const linal::Vector3D &MIN_BOX, const linal::Vector3D &MAX_BOX) const
      {
        const linal::Vector3D &mn( m_Box.GetMin()), &mx( m_Box.GetMax());
        double oob_dist_sq( 0.0);
        if( MAX_BOX.X() < mn.X())
        {
          oob_dist_sq += math::Sqr( MAX_BOX.X() - mn.X());
        }
        else if( MIN_BOX.X() > mx.X())
        {
          oob_dist_sq += math::Sqr( MIN_BOX.X() - mx.X());
        }
        if( MAX_BOX.Y() < mn.Y())
        {
          oob_dist_sq += math::Sqr( MAX_BOX.Y() - mn.Y());
        }
        else if( MIN_BOX.Y() > mx.Y())
        {
          oob_dist_sq += math::Sqr( MIN_BOX.Y() - mx.Y());
        }
        if( MAX_BOX.Z() < mn.Z())
        {
          oob_dist_sq += math::Sqr( MAX_BOX.Z() - mn.Z());
        }
        else if( MIN_BOX.Z() > mx.Z())
        {
          oob_dist_sq += math::Sqr( MIN_BOX.Z() - mx.Z());
        }
        return oob_dist_sq;
      }

      storage::Vector< storage::Triplet< SiPtr< const t_DataType>, SiPtr< const t_DataType>, double> >
        GetNeighborsMultiDimensional( const double &NEIGHBORHOOD) const
      {
        typedef storage::Triplet< SiPtr< const t_DataType>, SiPtr< const t_DataType>, double> t_Triplet;
        storage::Vector< t_Triplet> neighbors_ret;
        BCL_Assert
        (
          NEIGHBORHOOD <= m_Resolution,
          "GetNeighbors not supported for voxel grids with resolution smaller than neighborhood distance"
        );
        const double cutoff_sq( math::Sqr( NEIGHBORHOOD));

        if( m_CacheEdges)
        {
          auto itr_edges( m_Edges.Begin());
          for
          (
            auto itr_vox_obj( m_Assignments.Begin()), itr_vox_end( m_Assignments.End());
            itr_vox_obj != itr_vox_end;
            ++itr_vox_obj, ++itr_edges
          )
          {
            if( itr_vox_obj->IsEmpty())
            {
              continue;
            }

            auto itr_edges_begin( itr_edges->Begin()), itr_edges_end( itr_edges->End());
            // iterate over all neighbors within the voxel.
            // For high-density voxels (~32+ items) in 3D, it would be more efficient to do an
            // additional labeling of all items as to their sub-octant. Items in the voxel (on-the-fly) or as an auxiliary grid).
            // Any two objects with the same
            // octant number are automatically in the same neighborhood.
            for
            (
              auto itr_vox_obj_a( itr_vox_obj->Begin()), itr_vox_obj_a_end( itr_vox_obj->End());
              itr_vox_obj_a != itr_vox_obj_a_end;
              ++itr_vox_obj_a
            )
            {
              const linal::Vector3D &coord_a( *itr_vox_obj_a->m_Coordinates);
              for
              (
                auto itr_vox_obj_b( itr_vox_obj->Begin());
                itr_vox_obj_b != itr_vox_obj_a;
                ++itr_vox_obj_b
              )
              {
                const double sq_distance( linal::SquareDistance( coord_a, *itr_vox_obj_b->m_Coordinates));
                if( sq_distance < cutoff_sq)
                {
                  neighbors_ret.PushBack( t_Triplet( itr_vox_obj_a->m_Object, itr_vox_obj_b->m_Object, math::Sqrt( sq_distance)));
                }
              }

              for( auto itr_lists_b( itr_edges_begin); itr_lists_b != itr_edges_end; ++itr_lists_b)
              {
                for
                (
                  auto itr_vox_obj_b( ( *itr_lists_b)->Begin()), itr_vox_obj_b_end( ( *itr_lists_b)->End());
                  itr_vox_obj_b != itr_vox_obj_b_end;
                  ++itr_vox_obj_b
                )
                {
                  const double sq_distance( linal::SquareDistance( coord_a, *itr_vox_obj_b->m_Coordinates));
                  if( sq_distance < cutoff_sq)
                  {
                    neighbors_ret.PushBack( t_Triplet( itr_vox_obj_a->m_Object, itr_vox_obj_b->m_Object, math::Sqrt( sq_distance)));
                  }
                }
              }
            }
          }
        }
        else
        {
          size_t pos( 0);
          const int z_stride( m_NBinsX * m_NBinsY), nbx( m_NBinsX);
          auto itr_flag( m_BorderFlags.begin());
          for
          (
            auto itr_vox_obj( m_Assignments.Begin()), itr_vox_end( m_Assignments.End());
            itr_vox_obj != itr_vox_end;
            ++itr_vox_obj, ++pos, ++itr_flag
          )
          {
            if( itr_vox_obj->IsEmpty())
            {
              continue;
            }
            auto itr_vox_obj_a_end( itr_vox_obj->End());
            // iterate over all neighbors within the voxel.
            // For high-density voxels (~32+ items) in 3D, it would be more efficient to do an
            // additional labeling of all items as to their sub-octant. Items in the voxel (on-the-fly) or as an auxiliary grid).
            // Any two objects with the same
            // octant number are automatically in the same neighborhood.
            for( auto itr_vox_obj_a( itr_vox_obj->Begin()); itr_vox_obj_a != itr_vox_obj_a_end; ++itr_vox_obj_a)
            {
              const linal::Vector3D &coord_a( *itr_vox_obj_a->m_Coordinates);
              for
              (
                auto itr_vox_obj_b( itr_vox_obj->Begin());
                itr_vox_obj_b != itr_vox_obj_a;
                ++itr_vox_obj_b
              )
              {
                const double sq_distance( linal::SquareDistance( coord_a, *itr_vox_obj_b->m_Coordinates));
                if( sq_distance < cutoff_sq)
                {
                  neighbors_ret.PushBack( t_Triplet( itr_vox_obj_a->m_Object, itr_vox_obj_b->m_Object, math::Sqrt( sq_distance)));
                }
              }
            }

            const int flag( *itr_flag);
            const int max_z( flag & e_ZMax ? 1 : 2);
            const int min_y( flag & e_YMin ? 0 : -1), max_y( flag & e_YMax ? 1 : 2);
            const int min_x( flag & e_XMin ? 0 : -1), max_x( flag & e_XMax ? 1 : 2);
            for
            (
              int z( 0), zy_off_init( min_y * nbx), zy_off_end( max_y * nbx);
              z < max_z;
              ++z, zy_off_init += z_stride, zy_off_end += z_stride
            )
            {
              if( zy_off_end < 0)
              {
                continue;
              }
              for( int zy_off( zy_off_init); zy_off < zy_off_end; zy_off += nbx)
              {
                if( zy_off < 0)
                {
                  continue;
                }
                for( int zyx_off( zy_off + min_x), zyx_off_mx( zy_off + max_x); zyx_off < zyx_off_mx; ++zyx_off)
                {
                  if( zyx_off > 0)
                  {
                    const storage::Vector< ObjReference> &adjacent_vox( m_Assignments( pos + zyx_off));
                    if( adjacent_vox.IsEmpty())
                    {
                      continue;
                    }
                    auto itr_vox_obj_b_beg( adjacent_vox.Begin()), itr_vox_obj_b_end( adjacent_vox.End());
                    for( auto itr_vox_obj_a( itr_vox_obj->Begin()); itr_vox_obj_a != itr_vox_obj_a_end; ++itr_vox_obj_a)
                    {
                      const linal::Vector3D &coord_a( *itr_vox_obj_a->m_Coordinates);
                      for( auto itr_vox_obj_b( itr_vox_obj_b_beg); itr_vox_obj_b != itr_vox_obj_b_end; ++itr_vox_obj_b)
                      {
                        const double sq_distance( linal::SquareDistance( coord_a, *itr_vox_obj_b->m_Coordinates));
                        if( sq_distance < cutoff_sq)
                        {
                          neighbors_ret.PushBack( t_Triplet( itr_vox_obj_a->m_Object, itr_vox_obj_b->m_Object, math::Sqrt( sq_distance)));
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
        return neighbors_ret;
      }

      storage::Vector< storage::Pair< SiPtr< const t_DataType>, double> >
        GetNeighborsMultiDimensional( const t_DataType &OBJ, const double &NEIGHBORHOOD) const
      {
        typedef storage::Pair< SiPtr< const t_DataType>, double> t_Pair;
        storage::Vector< t_Pair> neighbors_ret;
        BCL_Assert
        (
          NEIGHBORHOOD <= m_Resolution,
          "GetNeighbors not supported for voxel grids with resolution smaller than neighborhood distance"
        );
        const double cutoff_sq( math::Sqr( NEIGHBORHOOD));
        SiPtr< const linal::Vector3D> coord_si( ExtractPosition( OBJ));
        if( !coord_si.IsDefined())
        {
          return neighbors_ret;
        }
        const linal::Vector3D &coord( *coord_si);
        if( GetSqDistanceOutsideBoundingBox( coord) >= cutoff_sq)
        {
          return neighbors_ret;
        }

        const size_t pos( GetIndex( coord));
        auto itr_vox_obj( m_Assignments.Begin() + pos);
        // iterate over all neighbors within the voxel
        for
        (
          auto itr_vox_obj_a( itr_vox_obj->Begin()), itr_vox_obj_vec_end( itr_vox_obj->End());
          itr_vox_obj_a != itr_vox_obj_vec_end;
          ++itr_vox_obj_a
        )
        {
          const double sq_distance( linal::SquareDistance( coord, *itr_vox_obj_a->m_Coordinates));
          if( sq_distance < cutoff_sq && !IsSameItem( *itr_vox_obj_a->m_Object, OBJ))
          {
            neighbors_ret.PushBack( t_Pair( itr_vox_obj_a->m_Object, math::Sqrt( sq_distance)));
          }
        }

        if( m_CacheEdges)
        {
          for
          (
            auto itr_lists_b( m_EdgesComplete( pos).Begin()), itr_lists_b_end( m_EdgesComplete( pos).End());
            itr_lists_b != itr_lists_b_end;
            ++itr_lists_b
          )
          {
            for
            (
              auto itr_vox_obj_b( ( *itr_lists_b)->Begin()), itr_vox_obj_b_end( ( *itr_lists_b)->End());
              itr_vox_obj_b != itr_vox_obj_b_end;
              ++itr_vox_obj_b
            )
            {
              const double sq_distance( linal::SquareDistance( coord, *itr_vox_obj_b->m_Coordinates));
              if( sq_distance < cutoff_sq)
              {
                neighbors_ret.PushBack( t_Pair( itr_vox_obj_b->m_Object, math::Sqrt( sq_distance)));
              }
            }
          }
        }
        else
        {
          // determine extents in each direction
          const int flag( m_BorderFlags[ pos]);
          const int min_z( flag & e_ZMin ? 0 : -1), max_z( flag & e_ZMax ? 1 : 2);
          const int min_y( flag & e_YMin ? 0 : -1), max_y( flag & e_YMax ? 1 : 2);
          const int min_x( flag & e_XMin ? 0 : -1), max_x( flag & e_XMax ? 1 : 2);
          const int z_stride( m_NBinsX * m_NBinsY);

          for
          (
            int z( min_z),
                zy_off_init( min_z * z_stride + min_y * m_NBinsX),
                zy_off_end( min_z * z_stride + max_y * m_NBinsX);
            z < max_z;
            ++z, zy_off_init += z_stride, zy_off_end += z_stride
          )
          {
            for( int zy_off( zy_off_init); zy_off < zy_off_end; zy_off += m_NBinsX)
            {
              for( int zyx_off( zy_off + min_x), zyx_off_mx( zy_off + max_x); zyx_off < zyx_off_mx; ++zyx_off)
              {
                if( zyx_off)
                {
                  const storage::Vector< ObjReference> &adjacent_vox( m_Assignments( pos + zyx_off));
                  for
                  (
                    auto itr_vox_obj_b( adjacent_vox.Begin()), itr_vox_obj_b_end( adjacent_vox.End());
                    itr_vox_obj_b != itr_vox_obj_b_end;
                    ++itr_vox_obj_b
                  )
                  {
                    const double sq_distance( linal::SquareDistance( coord, *itr_vox_obj_b->m_Coordinates));
                    if( sq_distance < cutoff_sq)
                    {
                      neighbors_ret.PushBack( t_Pair( itr_vox_obj_b->m_Object, math::Sqrt( sq_distance)));
                    }
                  }
                }
              }
            }
          }
        }
        return neighbors_ret;
      }

      storage::Vector< storage::Pair< SiPtr< const t_DataType>, double> >
        GetNeighborsMultiDimensional( const linal::Vector3D &COORD, const double &NEIGHBORHOOD) const
      {
        typedef storage::Pair< SiPtr< const t_DataType>, double> t_Pair;
        storage::Vector< t_Pair> neighbors_ret;
        BCL_Assert
        (
          NEIGHBORHOOD <= m_Resolution,
          "GetNeighbors not supported for voxel grids with resolution smaller than neighborhood distance"
        );
        const double cutoff_sq( math::Sqr( NEIGHBORHOOD));
        const linal::Vector3D &coord( COORD);
        if( GetSqDistanceOutsideBoundingBox( coord) >= cutoff_sq)
        {
          return neighbors_ret;
        }

        const size_t pos( GetIndex( coord));
        auto itr_vox_obj( m_Assignments.Begin() + pos);
        // iterate over all neighbors within the voxel
        for
        (
          auto itr_vox_obj_a( itr_vox_obj->Begin()), itr_vox_obj_vec_end( itr_vox_obj->End());
          itr_vox_obj_a != itr_vox_obj_vec_end;
          ++itr_vox_obj_a
        )
        {
          const double sq_distance( linal::SquareDistance( coord, *itr_vox_obj_a->m_Coordinates));
          if( sq_distance < cutoff_sq)
          {
            neighbors_ret.PushBack( t_Pair( itr_vox_obj_a->m_Object, math::Sqrt( sq_distance)));
          }
        }

        if( m_CacheEdges)
        {
          for
          (
            auto itr_lists_b( m_EdgesComplete( pos).Begin()), itr_lists_b_end( m_EdgesComplete( pos).End());
            itr_lists_b != itr_lists_b_end;
            ++itr_lists_b
          )
          {
            for
            (
              auto itr_vox_obj_b( ( *itr_lists_b)->Begin()), itr_vox_obj_b_end( ( *itr_lists_b)->End());
              itr_vox_obj_b != itr_vox_obj_b_end;
              ++itr_vox_obj_b
            )
            {
              const double sq_distance( linal::SquareDistance( coord, *itr_vox_obj_b->m_Coordinates));
              if( sq_distance < cutoff_sq)
              {
                neighbors_ret.PushBack( t_Pair( itr_vox_obj_b->m_Object, math::Sqrt( sq_distance)));
              }
            }
          }
        }
        else
        {
          // determine extents in each direction
          const int flag( m_BorderFlags[ pos]);
          const int min_z( flag & e_ZMin ? 0 : -1), max_z( flag & e_ZMax ? 1 : 2);
          const int min_y( flag & e_YMin ? 0 : -1), max_y( flag & e_YMax ? 1 : 2);
          const int min_x( flag & e_XMin ? 0 : -1), max_x( flag & e_XMax ? 1 : 2);
          const int z_stride( m_NBinsX * m_NBinsY);

          for
          (
            int z( min_z),
                zy_off_init( min_z * z_stride + min_y * m_NBinsX),
                zy_off_end( min_z * z_stride + max_y * m_NBinsX);
            z < max_z;
            ++z, zy_off_init += z_stride, zy_off_end += z_stride
          )
          {
            for( int zy_off( zy_off_init); zy_off < zy_off_end; zy_off += m_NBinsX)
            {
              for( int zyx_off( zy_off + min_x), zyx_off_mx( zy_off + max_x); zyx_off < zyx_off_mx; ++zyx_off)
              {
                if( zyx_off)
                {
                  const storage::Vector< ObjReference> &adjacent_vox( m_Assignments( pos + zyx_off));
                  for
                  (
                    auto itr_vox_obj_b( adjacent_vox.Begin()), itr_vox_obj_b_end( adjacent_vox.End());
                    itr_vox_obj_b != itr_vox_obj_b_end;
                    ++itr_vox_obj_b
                  )
                  {
                    const double sq_distance( linal::SquareDistance( coord, *itr_vox_obj_b->m_Coordinates));
                    if( sq_distance < cutoff_sq)
                    {
                      neighbors_ret.PushBack( t_Pair( itr_vox_obj_b->m_Object, math::Sqrt( sq_distance)));
                    }
                  }
                }
              }
            }
          }
        }
        return neighbors_ret;
      }

      storage::Vector< storage::Triplet< SiPtr< const t_DataType>, SiPtr< const t_DataType>, double> >
        GetNeighborsMultiDimensional( const VoxelGrid< t_DataType> &GRID, const double &NEIGHBORHOOD) const
      {
        typedef storage::Triplet< SiPtr< const t_DataType>, SiPtr< const t_DataType>, double> t_Triplet;
        storage::Vector< t_Triplet> neighbors_ret;
        BCL_Assert
        (
          NEIGHBORHOOD <= m_Resolution,
          "GetNeighbors not supported for voxel grids with resolution smaller than neighborhood distance"
        );
        const double cutoff_sq( math::Sqr( NEIGHBORHOOD));
        // test whether it is even possible that any points on the grids are within the cutoff distance of one another
        if( GetSqDistanceOutsideBoundingBox( GRID.m_Box.GetMin(), GRID.m_Box.GetMax()) > cutoff_sq)
        {
          return neighbors_ret;
        }

         for
        (
          auto itr_other_grid( GRID.m_Assignments.Begin()), itr_other_grid_end( GRID.m_Assignments.End());
          itr_other_grid != itr_other_grid_end;
          ++itr_other_grid
        )
        {
          for
          (
            auto itr_obj_a( itr_other_grid->Begin()), itr_obj_a_end( itr_other_grid->End());
            itr_obj_a != itr_obj_a_end;
            ++itr_obj_a
          )
          {
            const linal::Vector3D &coord( *itr_obj_a->m_Coordinates);
            if( GetSqDistanceOutsideBoundingBox( coord) >= cutoff_sq)
            {
              continue;
            }

            const SiPtr< const t_DataType> obj_a( itr_obj_a->m_Object);
            const size_t pos( GetIndex( coord));

            auto itr_vox_obj( m_Assignments.Begin() + pos);
            // iterate over all neighbors within the voxel
            for
            (
              auto itr_vox_obj_a( itr_vox_obj->Begin()), itr_vox_obj_vec_end( itr_vox_obj->End());
              itr_vox_obj_a != itr_vox_obj_vec_end;
              ++itr_vox_obj_a
            )
            {
              const double sq_distance( linal::SquareDistance( coord, *itr_vox_obj_a->m_Coordinates));
              if( sq_distance < cutoff_sq && !IsSameItem( *itr_vox_obj_a->m_Object, *obj_a))
              {
                neighbors_ret.PushBack( t_Triplet( itr_vox_obj_a->m_Object, obj_a, math::Sqrt( sq_distance)));
              }
            }

            if( m_CacheEdges)
            {
              for
              (
                auto itr_lists_b( m_EdgesComplete( pos).Begin()), itr_lists_b_end( m_EdgesComplete( pos).End());
                itr_lists_b != itr_lists_b_end;
                ++itr_lists_b
              )
              {
                for
                (
                  auto itr_vox_obj_b( ( *itr_lists_b)->Begin()), itr_vox_obj_b_end( ( *itr_lists_b)->End());
                  itr_vox_obj_b != itr_vox_obj_b_end;
                  ++itr_vox_obj_b
                )
                {
                  const double sq_distance( linal::SquareDistance( coord, *itr_vox_obj_b->m_Coordinates));
                  if( sq_distance < cutoff_sq)
                  {
                    neighbors_ret.PushBack( t_Triplet( itr_vox_obj_b->m_Object, obj_a, math::Sqrt( sq_distance)));
                  }
                }
              }
            }
            else
            {
              // determine extents in each direction
              const int flag( m_BorderFlags[ pos]);
              const int min_z( flag & e_ZMin ? 0 : -1), max_z( flag & e_ZMax ? 1 : 2);
              const int min_y( flag & e_YMin ? 0 : -1), max_y( flag & e_YMax ? 1 : 2);
              const int min_x( flag & e_XMin ? 0 : -1), max_x( flag & e_XMax ? 1 : 2);
              const int z_stride( m_NBinsX * m_NBinsY);

              for
              (
                int z( min_z),
                    zy_off_init( min_z * z_stride + min_y * m_NBinsX),
                    zy_off_end( min_z * z_stride + max_y * m_NBinsX);
                z < max_z;
                ++z, zy_off_init += z_stride, zy_off_end += z_stride
              )
              {
                for( int zy_off( zy_off_init); zy_off < zy_off_end; zy_off += m_NBinsX)
                {
                  for( int zyx_off( zy_off + min_x), zyx_off_mx( zy_off + max_x); zyx_off < zyx_off_mx; ++zyx_off)
                  {
                    if( zyx_off)
                    {
                      const storage::Vector< ObjReference> &adjacent_vox( m_Assignments( pos + zyx_off));
                      for
                      (
                        auto itr_vox_obj_b( adjacent_vox.Begin()), itr_vox_obj_b_end( adjacent_vox.End());
                        itr_vox_obj_b != itr_vox_obj_b_end;
                        ++itr_vox_obj_b
                      )
                      {
                        const double sq_distance( linal::SquareDistance( coord, *itr_vox_obj_b->m_Coordinates));
                        if( sq_distance < cutoff_sq)
                        {
                          neighbors_ret.PushBack( t_Triplet( itr_vox_obj_b->m_Object, obj_a, math::Sqrt( sq_distance)));
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
        return neighbors_ret;
      }

      storage::Vector< storage::Triplet< SiPtr< const t_DataType>, SiPtr< const t_DataType>, double> >
        GetNeighbors1D( const double &NEIGHBORHOOD) const
      {
        typedef storage::Triplet< SiPtr< const t_DataType>, SiPtr< const t_DataType>, double> t_Triplet;
        storage::Vector< t_Triplet> neighbors_ret;

        const double cutoff_sq( math::Sqr( NEIGHBORHOOD));
        size_t index_z_end( std::ceil( NEIGHBORHOOD / m_MinResolution) + 1);
        for
        (
          auto itr_zvec( m_Assignments.Begin()), itr_zvec_end( m_Assignments.End());
          itr_zvec != itr_zvec_end;
          ++itr_zvec, ++index_z_end
        )
        {
          auto itr_zvec_b_end( std::min( m_Assignments.Begin() + index_z_end, itr_zvec_end));
          for
          (
            auto itr_zobj( itr_zvec->Begin()), itr_zobj_end( itr_zvec->End());
            itr_zobj != itr_zobj_end;
            ++itr_zobj
          )
          {
            // for each object
            SiPtr< const t_DataType> p_obj( itr_zobj->m_Object);
            const linal::Vector3D &coord_ref( *itr_zobj->m_Coordinates);
            auto itr_zvec_b( itr_zvec);
            for( auto itr_short( itr_zvec_b->Begin()); itr_short != itr_zobj; ++itr_short)
            {
              const double dist_sq( linal::SquareDistance( coord_ref, *itr_short->m_Coordinates));
              if( dist_sq < cutoff_sq)
              {
                neighbors_ret.PushBack( t_Triplet( p_obj, itr_short->m_Object, math::Sqrt( dist_sq)));
              }
            }
            for( ++itr_zvec_b; itr_zvec_b != itr_zvec_b_end; ++itr_zvec_b)
            {
              for
              (
                auto itr_zobj_b( itr_zvec_b->Begin()), itr_zobj_b_end( itr_zvec_b->End());
                itr_zobj_b != itr_zobj_b_end;
                ++itr_zobj_b
              )
              {
                const double dist_sq( linal::SquareDistance( coord_ref, *itr_zobj_b->m_Coordinates));
                if( dist_sq < cutoff_sq)
                {
                  neighbors_ret.PushBack( t_Triplet( p_obj, itr_zobj_b->m_Object, math::Sqrt( dist_sq)));
                }
              }
            }
          }
        }
        return neighbors_ret;
      }

      storage::Vector< storage::Pair< SiPtr< const t_DataType>, double> >
        GetNeighbors1D( const t_DataType &OBJECT, const double &NEIGHBORHOOD) const
      {
        typedef storage::Pair< SiPtr< const t_DataType>, double> t_Pair;
        storage::Vector< t_Pair> neighbors_ret;

        SiPtr< const linal::Vector3D> coord_si( ExtractPosition( OBJECT));
        if( !coord_si.IsDefined())
        {
          return neighbors_ret;
        }
        const linal::Vector3D &coord( *coord_si);
        const double cutoff_sq( math::Sqr( NEIGHBORHOOD));
        if( GetSqDistanceOutsideBoundingBox( coord) >= cutoff_sq)
        {
          return neighbors_ret;
        }

        const size_t n_bins_1d( std::ceil( NEIGHBORHOOD / m_MinResolution));
        const size_t pos( GetIndex( coord));
        auto itr_zvec_b( m_Assignments.Begin() + pos), itr_zvec_b_base( m_Assignments.Begin() + pos);
        auto itr_zvec_b_end( std::min( m_Assignments.Begin() + pos + n_bins_1d + 1, m_Assignments.End()));
        for
        (
          auto itr_zvec_a( itr_zvec_b->Begin()), itr_zvec_a_end( itr_zvec_b->End());
          itr_zvec_a != itr_zvec_a_end;
          ++itr_zvec_a
        )
        {
          const double dist_sq( linal::SquareDistance( coord, *itr_zvec_a->m_Coordinates));
          if( dist_sq < cutoff_sq && !IsSameItem( *itr_zvec_a->m_Object, OBJECT))
          {
            neighbors_ret.PushBack( t_Pair( itr_zvec_a->m_Object, math::Sqrt( dist_sq)));
          }
        }
        for
        (
          itr_zvec_b = std::max( m_Assignments.Begin() - n_bins_1d + pos, m_Assignments.Begin());
          itr_zvec_b != itr_zvec_b_end;
          ++itr_zvec_b
        )
        {
          if( itr_zvec_b == itr_zvec_b_base)
          {
            continue;
          }
          for
          (
            auto itr_zobj_b( itr_zvec_b->Begin()), itr_zobj_b_end( itr_zvec_b->End());
            itr_zobj_b != itr_zobj_b_end;
            ++itr_zobj_b
          )
          {
            const double dist_sq( linal::SquareDistance( coord, *itr_zobj_b->m_Coordinates));
            if( dist_sq < cutoff_sq)
            {
              neighbors_ret.PushBack( t_Pair( itr_zobj_b->m_Object, math::Sqrt( dist_sq)));
            }
          }
        }
        return neighbors_ret;
      }

      storage::Vector< storage::Pair< SiPtr< const t_DataType>, double> >
        GetNeighbors1D( const linal::Vector3D &COORD, const double &NEIGHBORHOOD) const
      {
        typedef storage::Pair< SiPtr< const t_DataType>, double> t_Pair;
        storage::Vector< t_Pair> neighbors_ret;
        const linal::Vector3D &coord( COORD);
        const double cutoff_sq( math::Sqr( NEIGHBORHOOD));
        if( GetSqDistanceOutsideBoundingBox( coord) >= cutoff_sq)
        {
          return neighbors_ret;
        }

        const size_t n_bins_1d( std::ceil( NEIGHBORHOOD / m_MinResolution));
        const size_t pos( GetIndex( coord));
        auto itr_zvec_b( m_Assignments.Begin() + pos), itr_zvec_b_base( m_Assignments.Begin() + pos);
        auto itr_zvec_b_end( std::min( m_Assignments.Begin() + pos + n_bins_1d + 1, m_Assignments.End()));
        for
        (
          auto itr_zvec_a( itr_zvec_b->Begin()), itr_zvec_a_end( itr_zvec_b->End());
          itr_zvec_a != itr_zvec_a_end;
          ++itr_zvec_a
        )
        {
          const double dist_sq( linal::SquareDistance( coord, *itr_zvec_a->m_Coordinates));
          if( dist_sq < cutoff_sq)
          {
            neighbors_ret.PushBack( t_Pair( itr_zvec_a->m_Object, math::Sqrt( dist_sq)));
          }
        }
        for
        (
          itr_zvec_b = std::max( m_Assignments.Begin() - n_bins_1d + pos, m_Assignments.Begin());
          itr_zvec_b != itr_zvec_b_end;
          ++itr_zvec_b
        )
        {
          if( itr_zvec_b == itr_zvec_b_base)
          {
            continue;
          }
          for
          (
            auto itr_zobj_b( itr_zvec_b->Begin()), itr_zobj_b_end( itr_zvec_b->End());
            itr_zobj_b != itr_zobj_b_end;
            ++itr_zobj_b
          )
          {
            const double dist_sq( linal::SquareDistance( coord, *itr_zobj_b->m_Coordinates));
            if( dist_sq < cutoff_sq)
            {
              neighbors_ret.PushBack( t_Pair( itr_zobj_b->m_Object, math::Sqrt( dist_sq)));
            }
          }
        }
        return neighbors_ret;
      }

      storage::Vector< storage::Triplet< SiPtr< const t_DataType>, SiPtr< const t_DataType>, double> >
        GetNeighbors1D( const VoxelGrid< t_DataType> &GRID, const double &NEIGHBORHOOD) const
      {
        typedef storage::Triplet< SiPtr< const t_DataType>, SiPtr< const t_DataType>, double> t_Triplet;
        storage::Vector< t_Triplet> neighbors_ret;

        const double cutoff_sq( math::Sqr( NEIGHBORHOOD));
        // test whether it is even possible that any points on the grids are within the cutoff distance of one another
        if( GetSqDistanceOutsideBoundingBox( GRID.m_Box.GetMin(), GRID.m_Box.GetMax()) > cutoff_sq)
        {
          return neighbors_ret;
        }

        const size_t n_bins_1d( std::ceil( NEIGHBORHOOD / m_MinResolution));

        for
        (
          auto itr_zvec( GRID.m_Assignments.Begin()), itr_zvec_end( GRID.m_Assignments.End());
          itr_zvec != itr_zvec_end;
          ++itr_zvec
        )
        {
          for
          (
            auto itr_obj_a( itr_zvec->Begin()), itr_obj_a_end( itr_zvec->End());
            itr_obj_a != itr_obj_a_end;
            ++itr_obj_a
          )
          {
            const linal::Vector3D &coord( *itr_obj_a->m_Coordinates);
            if( GetSqDistanceOutsideBoundingBox( coord) >= cutoff_sq)
            {
              continue;
            }

            const SiPtr< const t_DataType> obj_a( itr_obj_a->m_Object);
            const size_t pos( GetIndex( coord));
            auto itr_zvec_b( m_Assignments.Begin() + pos), itr_zvec_b_base( m_Assignments.Begin() + pos);
            auto itr_zvec_b_end( std::min( m_Assignments.Begin() + pos + n_bins_1d + 1, m_Assignments.End()));
            for
            (
              auto itr_zvec_a( itr_zvec_b->Begin()), itr_zvec_a_end( itr_zvec_b->End());
              itr_zvec_a != itr_zvec_a_end;
              ++itr_zvec_a
            )
            {
              const double dist_sq( linal::SquareDistance( coord, *itr_zvec_a->m_Coordinates));
              if( dist_sq < cutoff_sq && !IsSameItem( *itr_zvec_a->m_Object, *obj_a))
              {
                neighbors_ret.PushBack( t_Triplet( itr_zvec_a->m_Object, obj_a, math::Sqrt( dist_sq)));
              }
            }
            for
            (
              itr_zvec_b = std::max( m_Assignments.Begin() - n_bins_1d + pos, m_Assignments.Begin());
              itr_zvec_b != itr_zvec_b_end;
              ++itr_zvec_b
            )
            {
              if( itr_zvec_b == itr_zvec_b_base)
              {
                continue;
              }
              for
              (
                auto itr_zobj_b( itr_zvec_b->Begin()), itr_zobj_b_end( itr_zvec_b->End());
                itr_zobj_b != itr_zobj_b_end;
                ++itr_zobj_b
              )
              {
                const double dist_sq( linal::SquareDistance( coord, *itr_zobj_b->m_Coordinates));
                if( dist_sq < cutoff_sq)
                {
                  neighbors_ret.PushBack( t_Triplet( itr_zobj_b->m_Object, obj_a, math::Sqrt( dist_sq)));
                }
              }
            }
          }
        }
        return neighbors_ret;
      }

    }; //class VoxelGrid
  } // namespace util
} // namespace bcl

#endif //BCL_UTIL_VOXEL_GRID_H_

#pragma once

#include <memory>
#include <utility>
#include "util/mempool.h"

namespace seal
{
    /**
    SEAL uses memory pools for improved performance due to the large number of memory allocations needed
    by the homomorphic encryption operations, and the underlying polynomial arithmetic. The library 
    automatically creates a shared global memory pool, that is by default used by all instances of the 
    computation-heavy classes such as Encryptor, Evaluator, and PolyCRTBuilder. However, sometimes the 
    user might want to use local memory pools with some of these classes. For example, in heavily 
    multi-threaded applications the global memory pool might become clogged due to concurrent allocations. 
    Instead, the user might want to create a separate---say Evaluator---object for each thread, and have it
    use a thread-local memory pool. The MemoryPoolHandle class provides the functionality for doing this.

    For example, the user can create a MemoryPoolHandle that points to a new local memory pool by calling
    the static function acquire_new(). The MemoryPoolHandle it returns (or copies of it) can now be passed 
    on as an argument to the constructors of one or more classes (such as Encryptor, Evaluator, and PolyCRTBuilder). 

    Internally, a MemoryPoolHandle simply wraps a shared pointer to a util::MemoryPool object. This way
    the local memory pool will be automatically destroyed, and the memory released, as soon as no existing 
    handles point to it. Since the global memory pool is a static object, it will always have a positive 
    reference count, and thus will not be destroyed until the program terminates.
    */
    class MemoryPoolHandle
    {
    public:
        /**
        Creates a new MemoryPoolHandle pointing to the global memory pool.
        */
        MemoryPoolHandle() : pool_(util::MemoryPoolMT::default_pool())
        {
        }

        /**
        Overwrites the MemoryPoolHandle instance with the specified instance, so that the current instance 
        will point to the same underlying memory pool as the assigned instance.

        @param[in] assign The MemoryPoolHandle instance that should be assigned to the current instance
        */
        inline MemoryPoolHandle &operator =(const MemoryPoolHandle &assign)
        {
            pool_ = assign.pool_;
            return *this;
        }

        /**
        Overwrites the MemoryPoolHandle instance with the specified instance, so that the current instance 
        will point to the same underlying memory pool as the assigned instance.

        @param[in] assign The MemoryPoolHandle instance that should be assigned to the current instance
        */
        inline MemoryPoolHandle &operator =(MemoryPoolHandle &&assign) noexcept
        {
            pool_ = std::move(assign.pool_);
            return *this;
        }

        /**
        Creates a copy of a MemoryPoolHandle.

        @param[in] copy The MemoryPoolHandle to copy from
        */
        MemoryPoolHandle(const MemoryPoolHandle &copy)
        {
            operator =(copy);
        }

        /**
        Creates a new MemoryPoolHandle by moving an old one.

        @param[in] source The MemoryPoolHandle to move from
        */
        MemoryPoolHandle(MemoryPoolHandle &&source) noexcept
        {
            operator =(std::move(source));
        }

        /**
        Returns a MemoryPoolHandle pointing to the global memory pool.
        */
        inline static MemoryPoolHandle acquire_global()
        {
            return MemoryPoolHandle();
        }

        /**
        Returns a MemoryPoolHandle pointing to a new thread-safe memory pool. Optionally, the user
        can set the parameter thread_safe to false to create a thread-unsafe memory pool, which
        might be faster in certain use-cases where only single-threaded execution is required.

        @param[in] thread_safe Whether to create a thread-safe (default) or thread-unsafe memory pool
        */
        inline static MemoryPoolHandle acquire_new(bool thread_safe = true)
        {
            if (thread_safe)
            {
                return MemoryPoolHandle(std::make_shared<util::MemoryPoolMT>());
            }
            return MemoryPoolHandle(std::make_shared<util::MemoryPoolST>());
        }

        /**
        Returns a reference to the internal util::MemoryPool that the MemoryPoolHandle points to.
        */
        inline operator util::MemoryPool &() const
        {
            return *pool_.get();
        }

        /**
        Returns the number of pools (different size allocations) the internal util::MemoryPool has made. 
        */        
        inline std::uint64_t pool_count() const
        {
            return pool_->pool_count();
        }

        /**
        Returns the total amount of memory the internal util::MemoryPool has allocated.
        */
        inline std::uint64_t alloc_uint64_count() const
        {
            return pool_->alloc_uint64_count();
        }

        /**
        Returns the total amount of memory the internal util::MemoryPool has allocated.
        */
        inline std::uint64_t alloc_byte_count() const
        {
            return pool_->alloc_byte_count();
        }

    private:
        MemoryPoolHandle(std::shared_ptr<util::MemoryPool> &&pool) noexcept : pool_(std::move(pool))
        {
        }

        std::shared_ptr<util::MemoryPool> pool_;
    };
}
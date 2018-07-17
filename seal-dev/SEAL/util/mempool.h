#pragma once

#include <cstdint>
#include <vector>
#include <stdexcept>
#include <memory>
#include "util/common.h"
#include "util/locks.h"

namespace seal
{
    namespace util
    {
        class MemoryPoolItem
        {
        public:
            MemoryPoolItem(std::uint64_t *pointer) : pointer_(pointer), next_(nullptr)
            {
            }

            std::uint64_t *pointer()
            {
                return pointer_;
            }

            const std::uint64_t *pointer() const
            {
                return pointer_;
            }

            MemoryPoolItem* &next()
            {
                return next_;
            }

            const MemoryPoolItem *next() const
            {
                return next_;
            }

        private:
            MemoryPoolItem(const MemoryPoolItem &copy) = delete;

            MemoryPoolItem &operator =(const MemoryPoolItem &assign) = delete;

            std::uint64_t *pointer_;

            MemoryPoolItem *next_;
        };

        class MemoryPoolHead
        {
        public:
            struct allocation
            {
                static const std::uint64_t first_alloc_count;

                static const double alloc_size_multiplier;

                // Size of an allocation (in uint64_count)
                std::uint64_t size;

                // Pointer to start of allocation
                std::uint64_t *ptr;

                // How much free space is left (in uint64_count)
                std::uint64_t free;

                // Pointer to current head of allocation
                std::uint64_t *head_ptr;
            };

            virtual ~MemoryPoolHead()
            {
            }

            virtual std::uint64_t uint64_count() const = 0;

            virtual std::uint64_t alloc_item_count() const = 0;

            virtual MemoryPoolItem *get() = 0;

            virtual void add(MemoryPoolItem *new_first) = 0;
        };

        class MemoryPoolHeadMT : public MemoryPoolHead
        {
        public:
            // Creates a new MemoryPoolHeadMT with allocation for one single item.
            MemoryPoolHeadMT(std::uint64_t uint64_count);

            ~MemoryPoolHeadMT() override;

            std::uint64_t uint64_count() const override
            {
                return uint64_count_;
            }

            // Returns the total number of items allocated
            std::uint64_t alloc_item_count() const override
            {
                return alloc_item_count_;
            }

            MemoryPoolItem *get() override;

            void add(MemoryPoolItem *new_first) override;

        private:
            MemoryPoolHeadMT(const MemoryPoolHeadMT &copy) = delete;

            MemoryPoolHeadMT &operator =(const MemoryPoolHeadMT &assign) = delete;

            mutable std::atomic<bool> locked_;

            volatile std::uint64_t uint64_count_;

            volatile std::uint64_t alloc_item_count_;

            std::vector<allocation> allocs_;

            MemoryPoolItem* volatile first_item_;
        };

        class MemoryPoolHeadST : public MemoryPoolHead
        {
        public:
            // Creates a new MemoryPoolHeadST with allocation for one single item.
            MemoryPoolHeadST(std::uint64_t uint64_count);

            ~MemoryPoolHeadST() override;

            std::uint64_t uint64_count() const override
            {
                return uint64_count_;
            }

            // Returns the total number of items allocated
            std::uint64_t alloc_item_count() const override
            {
                return alloc_item_count_;
            }

            MemoryPoolItem *get() override;

            inline void add(MemoryPoolItem* new_first) override
            {
                new_first->next() = first_item_;
                first_item_ = new_first;
            }

        private:
            MemoryPoolHeadST(const MemoryPoolHeadST &copy) = delete;

            MemoryPoolHeadST &operator =(const MemoryPoolHeadST &assign) = delete;

            std::uint64_t uint64_count_;

            std::uint64_t alloc_item_count_;

            std::vector<allocation> allocs_;

            MemoryPoolItem *first_item_;
        };

        class ConstPointer;

        class Pointer
        {
        public:
            friend class ConstPointer;

            Pointer() : pointer_(nullptr), head_(nullptr), item_(nullptr), alias_(false)
            {
            }

            Pointer(MemoryPoolHead *head) : pointer_(nullptr), head_(nullptr), item_(nullptr), alias_(false)
            {
#ifdef _DEBUG
                if (head == nullptr)
                {
                    throw std::invalid_argument("head");
                }
#endif
                head_ = head;
                item_ = head->get();
                pointer_ = item_->pointer();
            }

            Pointer(Pointer &&move) noexcept : pointer_(move.pointer_), head_(move.head_), item_(move.item_), alias_(move.alias_)
            {
                move.pointer_ = nullptr;
                move.head_ = nullptr;
                move.item_ = nullptr;
                move.alias_ = false;
            }

            std::uint64_t &operator[](int index)
            {
                return pointer_[index];
            }

            std::uint64_t operator[](int index) const
            {
                return pointer_[index];
            }

            Pointer &operator =(Pointer &&assign)
            {
                acquire(assign);
                return *this;
            }

            bool is_set() const
            {
                return pointer_ != nullptr;
            }

            std::uint64_t *get()
            {
                return pointer_;
            }

            const std::uint64_t *get() const
            {
                return pointer_;
            }

            void release()
            {
                if (head_ != nullptr)
                {
                    head_->add(item_);
                }
                else if (pointer_ != nullptr && !alias_)
                {
                    delete[] pointer_;
                }
                pointer_ = nullptr;
                head_ = nullptr;
                item_ = nullptr;
                alias_ = false;
            }

            void acquire(Pointer &other)
            {
                if (this == &other)
                {
                    return;
                }
                if (head_ != nullptr)
                {
                    head_->add(item_);
                }
                else if (pointer_ != nullptr && !alias_)
                {
                    delete[] pointer_;
                }
                pointer_ = other.pointer_;
                head_ = other.head_;
                item_ = other.item_;
                alias_ = other.alias_;
                other.pointer_ = nullptr;
                other.head_ = nullptr;
                other.item_ = nullptr;
                other.alias_ = false;
            }

            void swap_with(Pointer &other) noexcept
            {
                std::swap(pointer_, other.pointer_);
                std::swap(head_, other.head_);
                std::swap(item_, other.item_);
                std::swap(alias_, other.alias_);
            }

            ~Pointer()
            {
                release();
            }

            static Pointer Owning(std::uint64_t *pointer)
            {
                return Pointer(pointer, false);
            }

            static Pointer Aliasing(std::uint64_t *pointer)
            {
                return Pointer(pointer, true);
            }

        private:
            Pointer(std::uint64_t *pointer, bool alias) : pointer_(pointer), head_(nullptr), item_(nullptr), alias_(alias)
            {
            }

            Pointer(const Pointer &copy) = delete;

            Pointer &operator =(const Pointer &assign) = delete;

            std::uint64_t *pointer_;

            MemoryPoolHead *head_;

            MemoryPoolItem *item_;

            bool alias_;
        };

        class ConstPointer
        {
        public:
            ConstPointer() : pointer_(nullptr), head_(nullptr), item_(nullptr), alias_(false)
            {
            }

            ConstPointer(MemoryPoolHead *head) : pointer_(nullptr), head_(nullptr), item_(nullptr), alias_(false)
            {
#ifdef _DEBUG
                if (head == nullptr)
                {
                    throw std::invalid_argument("head");
                }
#endif
                pointer_ = item_->pointer();
                head_ = head;
                item_ = head->get();
            }

            ConstPointer(ConstPointer &&move) noexcept : pointer_(move.pointer_), head_(move.head_), item_(move.item_), alias_(move.alias_)
            {
                move.pointer_ = nullptr;
                move.head_ = nullptr;
                move.item_ = nullptr;
                move.alias_ = false;
            }

            ConstPointer(Pointer &&move) noexcept : pointer_(move.pointer_), head_(move.head_), item_(move.item_), alias_(move.alias_)
            {
                move.pointer_ = nullptr;
                move.head_ = nullptr;
                move.item_ = nullptr;
                move.alias_ = false;
            }

            ConstPointer &operator =(ConstPointer &&assign)
            {
                acquire(assign);
                return *this;
            }

            ConstPointer &operator =(Pointer &&assign)
            {
                acquire(assign);
                return *this;
            }

            std::uint64_t operator[](int index) const
            {
                return pointer_[index];
            }

            bool is_set() const
            {
                return pointer_ != nullptr;
            }

            const std::uint64_t *get() const
            {
                return pointer_;
            }

            void release()
            {
                if (head_ != nullptr)
                {
                    head_->add(item_);
                }
                else if (pointer_ != nullptr && !alias_)
                {
                    delete[] pointer_;
                }
                pointer_ = nullptr;
                head_ = nullptr;
                item_ = nullptr;
                alias_ = false;
            }

            void acquire(ConstPointer &other)
            {
                if (this == &other)
                {
                    return;
                }
                if (head_ != nullptr)
                {
                    head_->add(item_);
                }
                else if (pointer_ != nullptr && !alias_)
                {
                    delete[] pointer_;
                }                
                pointer_ = other.pointer_;
                head_ = other.head_;
                item_ = other.item_;
                alias_ = other.alias_;
                other.pointer_ = nullptr;
                other.head_ = nullptr;
                other.item_ = nullptr;
                other.alias_ = false;
            }

            void acquire(Pointer &other)
            {
                if (head_ != nullptr)
                {
                    head_->add(item_);
                }
                else if (pointer_ != nullptr && !alias_)
                {
                    delete[] pointer_;
                }                
                pointer_ = other.pointer_;
                head_ = other.head_;
                item_ = other.item_;
                alias_ = other.alias_;
                other.pointer_ = nullptr;
                other.head_ = nullptr;
                other.item_ = nullptr;
                other.alias_ = false;
            }

            void swap_with(ConstPointer &other) noexcept
            {
                std::swap(pointer_, other.pointer_);
                std::swap(head_, other.head_);
                std::swap(item_, other.item_);
                std::swap(alias_, other.alias_);
            }

            ~ConstPointer()
            {
                release();
            }

            static ConstPointer Owning(std::uint64_t *pointer)
            {
                return ConstPointer(pointer, false);
            }

            static ConstPointer Aliasing(const std::uint64_t *pointer)
            {
                return ConstPointer(const_cast<uint64_t*>(pointer), true);
            }

        private:
            ConstPointer(std::uint64_t *pointer, bool alias) : pointer_(pointer), head_(nullptr), item_(nullptr), alias_(alias)
            {
            }

            ConstPointer(const ConstPointer &copy) = delete;

            ConstPointer &operator =(const ConstPointer &assign) = delete;

            std::uint64_t *pointer_;

            MemoryPoolHead *head_;

            MemoryPoolItem *item_;

            bool alias_;
        };

        class MemoryPool
        {
        public:
            virtual ~MemoryPool()
            {
            }

            virtual Pointer get_for_byte_count(std::uint64_t byte_count) = 0;

            virtual Pointer get_for_uint64_count(std::uint64_t uint64_count) = 0;

            virtual std::uint64_t pool_count() const = 0;

            virtual std::uint64_t alloc_uint64_count() const = 0;

            virtual std::uint64_t alloc_byte_count() const = 0;
        };

        class MemoryPoolMT : public MemoryPool
        {
        public:
            MemoryPoolMT() 
            {
            }

            ~MemoryPoolMT();

            Pointer get_for_byte_count(std::uint64_t byte_count)
            {
                std::uint64_t uint64_count = (byte_count + bytes_per_uint64 - 1) / bytes_per_uint64;
                return get_for_uint64_count(uint64_count);
            }

            Pointer get_for_uint64_count(std::uint64_t uint64_count);

            std::uint64_t pool_count() const
            {
                ReaderLock lock = pools_locker_.acquire_read();
                return pools_.size();
            }

            std::uint64_t alloc_uint64_count() const;

            std::uint64_t alloc_byte_count() const
            {
                return alloc_uint64_count() * bytes_per_uint64;
            }

            static std::shared_ptr<MemoryPoolMT> default_pool()
            {
                static std::shared_ptr<MemoryPoolMT> default_pool_(std::make_shared<MemoryPoolMT>());

                return default_pool_;
            }

        private:
            MemoryPoolMT(const MemoryPoolMT &copy) = delete;

            MemoryPoolMT &operator =(const MemoryPoolMT &assign) = delete;

            mutable ReaderWriterLocker pools_locker_;

            std::vector<MemoryPoolHead*> pools_;
        };

        class MemoryPoolST : public MemoryPool
        {
        public:
            MemoryPoolST()
            {
            }

            ~MemoryPoolST();

            Pointer get_for_byte_count(std::uint64_t byte_count)
            {
                std::uint64_t uint64_count = (byte_count + bytes_per_uint64 - 1) / bytes_per_uint64;
                return get_for_uint64_count(uint64_count);
            }

            Pointer get_for_uint64_count(std::uint64_t uint64_count);

            std::uint64_t pool_count() const
            {
                return pools_.size();
            }

            std::uint64_t alloc_uint64_count() const;

            std::uint64_t alloc_byte_count() const
            {
                return alloc_uint64_count() * bytes_per_uint64;
            }

        private:
            MemoryPoolST(const MemoryPoolST &copy) = delete;

            MemoryPoolST &operator =(const MemoryPoolST &assign) = delete;

            std::vector<MemoryPoolHead*> pools_;
        };

        Pointer duplicate_if_needed(std::uint64_t *original, int uint64_count, bool condition, MemoryPool &pool);

        ConstPointer duplicate_if_needed(const std::uint64_t *original, int uint64_count, bool condition, MemoryPool &pool);
    }
}

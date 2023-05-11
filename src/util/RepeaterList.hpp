//
// Created by juli on 31.07.22.
//

#ifndef GENEPISEEKER_REPEATERLIST_HPP
#define GENEPISEEKER_REPEATERLIST_HPP

#include <vector>
#include <memory>
#include <stdexcept>

// #define REPEATER_LIST_PRINT

namespace epi {
    /*!
     * A vector like container where one can assign an amount (integer) for every item,
     * The items can be accessed by index as if all elements would be in the list the given number of times.
     * Example: [0,0,0,1,1,2,3,3] --> [{3 x 0}, {2 x 1}, {1 x 2}, {2 x 3}]
     * Can only be created from a vector of pairs containing the element and the amount. Elements can be deleted efficiently by index.
     * No new elements can be added. Item amounts can, however, be restored if one knows their id
     * @tparam T
     */
    template<class T>
    class RepeaterList {
        class TreeItem {
        public:
            static std::shared_ptr<TreeItem> create(typename std::vector<std::pair<T, std::size_t>>::iterator begin, typename std::vector<std::pair<T, std::size_t>>::iterator end) {
                size_t num_items = end - begin;
                if (num_items == 0) {
                    return std::make_shared<TreeItemEmpty>();
                } else if (num_items == 1) {
                    return std::make_shared<TreeItemLeaf>(*begin);
                } else {
                    // more than two objects --> split into two subparts
                    size_t mid = num_items / 2;
                    return std::make_shared<TreeItemNode>(create(begin, begin + mid), create(begin + mid, end));
                }
            }

            virtual size_t build_index() = 0;
            virtual size_t build_id() = 0;
            virtual T at(std::size_t pos) = 0;
            virtual size_t remove_at(std::size_t pos) = 0;
            virtual std::size_t remove_group_at(std::size_t pos) = 0;
            virtual std::size_t get_group_end(std::size_t pos) = 0;
            virtual std::size_t get_group_start(std::size_t pos) = 0;

            virtual void restore_group_item(std::size_t group_id) = 0;

#ifdef REPEATER_LIST_PRINT
            virtual void printInorder(std::ostream & os) = 0;
#endif

            virtual ~TreeItem() = default;
        };

        class TreeItemNode : public TreeItem {
            std::shared_ptr<TreeItem> left;
            std::shared_ptr<TreeItem> right;

            size_t separator = 0;
            size_t id_separator = 0;

        public:
            TreeItemNode(std::shared_ptr<TreeItem> _left, std::shared_ptr<TreeItem> _right) {
                left = _left;
                right = _right;
            }

            size_t build_index() override {
                separator = left->build_index();
                return right->build_index() + separator;
            }

            size_t build_id() override {
                id_separator = left->build_id();
                return right->build_id() + id_separator;
            }

            void restore_group_item(std::size_t group_id) override {
                if (group_id < id_separator) {
                    left->restore_group_item(group_id);
                    separator++;
                } else {
                    right->restore_group_item(group_id - id_separator);
                }
            }

            T at(std::size_t pos) override {
                if (pos < separator) return left->at(pos);
                return right->at(pos - separator);
            }

            size_t remove_at(std::size_t pos) override {
                if (pos < separator) {
                    auto id = left->remove_at(pos);
                    separator --;
                    return id;
                } else {
                    return right->remove_at(pos - separator) + id_separator;
                }
            }

            std::size_t remove_group_at(std::size_t pos) override {
                if (pos < separator) {
                    auto num_removed = left->remove_group_at(pos);
                    separator -= num_removed;
                    return num_removed;
                } else {
                    return right->remove_group_at(pos - separator);
                }
            }

            std::size_t get_group_end(std::size_t pos) override {
                if (pos < separator) {
                    return left->get_group_end(pos);
                } else {
                    return right->get_group_end(pos - separator) + separator;
                }
            }

            std::size_t get_group_start(std::size_t pos) override {
                if (pos < separator) {
                    return left->get_group_start(pos);
                } else {
                    return right->get_group_start(pos - separator) + separator;
                }
            }

#ifdef REPEATER_LIST_PRINT
            void printInorder(std::ostream &os) override {
                left->printInorder(os);
                os << ", ";
                right->printInorder(os);
            }
#endif
        };

        class TreeItemLeaf : public TreeItem {
            T data;
            size_t amount;
        public:
            explicit TreeItemLeaf(std::pair<T, std::size_t> item) {
                data = item.first;
                amount = item.second;
            }

            size_t build_index () override {
                return amount;
            }

            size_t build_id () override {
                return 1;
            }

            T at(std::size_t pos) override {
                return data;
            }

            size_t remove_at(std::size_t pos) override {
                amount--;
                return 0;
            }

            void restore_group_item(std::size_t group_id) override {
                amount++;
            }

            std::size_t remove_group_at(std::size_t pos) override {
                size_t prev = amount;
                amount = 0;
                return prev;
            }

            std::size_t get_group_end(std::size_t pos) override {
                return amount - 1;
            }

            std::size_t get_group_start(std::size_t pos) override {
                return 0;
            }

#ifdef REPEATER_LIST_PRINT
            void printInorder(std::ostream &os) override {
                os << "{ " << amount << " x '" << data << "' }";
            }
#endif
        };

        class TreeItemEmpty : public TreeItem {
        public:
            size_t build_index() override {
                return 0;
            }

            size_t build_id() override {
                return 0;
            }

            T at(std::size_t pos) override {
                return {};
            }

            size_t remove_at(std::size_t pos) override {
                return 0;
            }
            void restore_group_item(std::size_t group_id) override {}
            std::size_t remove_group_at(std::size_t pos) override {
                return 0;
            }

            std::size_t get_group_end(std::size_t pos) override {
                return 0;
            }

            std::size_t get_group_start(std::size_t pos) override {
                return 0;
            }

#ifdef REPEATER_LIST_PRINT
            void printInorder(std::ostream &os) override {}
#endif
        };


        std::shared_ptr<TreeItem> root;
        size_t _size;
        size_t _distinct_size;
    public:
        explicit RepeaterList(std::vector<std::pair<T, std::size_t>> input) {
            root = TreeItem::create(input.begin(), input.end());
            _size = root->build_index();
            _distinct_size = root->build_id();
        }
        T operator[](std::size_t pos) {
            if (pos >= _size) {
                throw std::logic_error("Tried to access invalid index " + std::to_string(pos) + " (size = " + std::to_string(_size) + ")");
            }
            return root->at(pos);
        }

        /**
         * erases one element at the given position and returns the group id (needed to restore this action later)
         * @param pos
         * @return group id (needed to restore the action)
         */
        std::size_t erase(std::size_t pos) {
            if (pos >= _size) {
                throw std::logic_error("Tried to erase invalid index " + std::to_string(pos) + " (size = " + std::to_string(_size) + ")");
            }

            auto id = root->remove_at(pos);
            _size --;
            return id;
        }

        std::size_t erase_group(std::size_t pos) {
            if (pos >= _size) {
                throw std::logic_error("Tried to erase invalid index " + std::to_string(pos) + " (size = " + std::to_string(_size) + ")");
            }

            size_t num_removed = root->remove_group_at(pos);
            _size -= num_removed;
            return num_removed;
        }

        size_t get_group_end(std::size_t pos) {
            if (pos >= _size) {
                throw std::logic_error("Tried to access invalid index " + std::to_string(pos) + " (size = " + std::to_string(_size) + ")");
            }

            return root->get_group_end(pos);
        }

        size_t get_group_start(std::size_t pos) {
            if (pos >= _size) {
                throw std::logic_error("Tried to access invalid index " + std::to_string(pos) + " (size = " + std::to_string(_size) + ")");
            }

            return root->get_group_start(pos);
        }

        void restore_item_of_group(std::size_t group_id) {
            if (group_id >= _distinct_size) {
                throw std::logic_error("Tried to restore an item of a group with invalid id " + std::to_string(group_id) + " (#groups = " + std::to_string(_distinct_size) + ")");
            }

            root->restore_group_item(group_id);
            _size++;
        }

        size_t size() {
            return _size;
        }

#ifdef REPEATER_LIST_PRINT
        friend std::ostream & operator<<(std::ostream & os, const RepeaterList<T>& li) {
            os << '[';
            li.root->printInorder(os);
            os << "]  (total: " << li._size << ", distinct: " << li._distinct_size << ')';
            return os;
        }
#endif
    };
} // epi

#endif //GENEPISEEKER_REPEATERLIST_HPP

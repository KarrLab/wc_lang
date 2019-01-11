
    def merge_models(self, primary_model, secondary_model):
        """ Merge two models by merging objects from secondary model into primary model

        Args:
            primary_model (:obj:`Model`): base for merged model
            secondary_model (:obj:`Model`): model to merge into primary model

        Raises:
            :obj:`ValueError`: if an unsupported merge operation is encountered
            :obj:`ValueError`: if two objects with the same primary attribute cannot be merged because objects of their class cannot be merged
        """
        shared_objs, join_objs, different_secondary_objs = self.get_common_uncommon_objs(primary_model, secondary_model)

        # merge objects that have the same primary attributes
        for secondary_obj, primary_obj in join_objs.items():
            merge_type = primary_obj.Meta.merge
            if merge_type == MergeType.primary:
                continue
            elif merge_type == MergeType.join:
                self.join_objs(primary_model, primary_obj, secondary_model, secondary_obj, shared_objs)
            elif merge_type == MergeType.append:
                raise ValueError('Cannot merge {} {} from {} and {}'.format(
                    primary_obj.Meta.verbose_name, primary_obj.get_primary_attribute(),
                    primary_model.id, secondary_model.id))
            else:
                raise ValueError('Unsupported merge operation: "{}"'.format(merge_type.name)
                                 )  # pragma: no cover # unreachable because all values of MergeType are enumerated

        for different_secondary_obj in different_secondary_objs:
            self.append_obj(different_secondary_obj, shared_objs)

    def get_common_uncommon_objs(self, primary_model, secondary_model):
        """ Get the intersection and difference of the related objects of the primary and secondary models

        Args:
            primary_model (:obj:`Model`): base for merged model
            secondary_model (:obj:`Model`): model to merge into primary model

        Returns:
            :obj:`dict`: dictionary that maps objects in the secondary model to objects in the primary
                model
            :obj:`dict`: dictionary that maps objects in the secondary model that should be joined with
                objects in the primary model
            :obj:`list`: list of objects in the secondary model that have different primary attributes (e.g. id)
                from the objects in the primary model

        Raises:
            :obj:`ValueError`: if instances of the same class with the same primary attribute cannot be merged
                because the objects of the class should be inherited from the primary model and there are multiple
                instances of the class
        """
        primary_objs = primary_model.get_related_objects()
        secondary_objs = secondary_model.get_related_objects()

        # identify common and different objects
        primary_objs_by_class = {}
        secondary_objs_by_class = {}
        shared_objs = {}
        join_objs = {}
        different_secondary_objs = []
        for model in get_models():
            primary_objs_by_class[model] = {}
            secondary_objs_by_class[model] = {}

        for primary_obj in primary_objs:
            serialized_val = primary_obj.serialize()
            primary_objs_by_class[primary_obj.__class__][serialized_val] = primary_obj

        for secondary_obj in secondary_objs:
            serialized_val = secondary_obj.serialize()
            secondary_objs_by_class[secondary_obj.__class__][serialized_val] = secondary_obj
            if secondary_obj.Meta.merge == MergeType.primary:
                if not primary_objs_by_class[primary_obj.__class__]:
                    shared_objs[secondary_obj] = None
                elif len(primary_objs_by_class[primary_obj.__class__]) == 1:
                    shared_objs[secondary_obj] = list(primary_objs_by_class[primary_obj.__class__].values())[0]
                else:
                    raise ValueError('Unable to merge {} "{}"'.format(secondary_obj.Meta.verbose_name.lower(), secondary_obj.serialize()))

            elif serialized_val in primary_objs_by_class[primary_obj.__class__]:
                shared_objs[secondary_obj] = primary_objs_by_class[primary_obj.__class__][serialized_val]
                join_objs[secondary_obj] = primary_objs_by_class[primary_obj.__class__][serialized_val]

            else:
                different_secondary_objs.append(secondary_obj)

        return (shared_objs, join_objs, different_secondary_objs)

    def join_objs(self, primary_model, primary_obj, secondary_model, secondary_obj, shared_objs):
        """ Join two objects from the primary and secondary models

        Args:
            primary_model (:obj:`Model`): base for merged model
            primary_obj (): object from primary model to merge with :obj:`secondary_obj`
            secondary_model (:obj:`Model`): model to merge into primary model
            secondary_obj (:obj:`obj_model.Model`): object from secondary model to merge with :obj:`primary_obj`
            shared_objects (:obj:`dict`): dictionary that maps objects in the secondary model to objects in the primary
                model that have the same primary attribute (e.g. id)

        Raises:
            :obj:`ValueError`: if the primary and secondary objects have differt values of a scalar attribute
            :obj:`ValueError`: if the merge operation is not supported
        """
        for attr in primary_obj.Meta.attributes.values():
            primary_val = getattr(primary_obj, attr.name)
            secondary_val = getattr(secondary_obj, attr.name)

            if isinstance(attr, obj_model.LiteralAttribute):
                if primary_val != secondary_val:
                    raise ValueError('{} {} of {} and {} must have the same value of {}'.format(
                        primary_obj.Meta.verbose_name, primary_obj.get_primary_attribute(),
                        primary_model.id, secondary_model.id, attr.name))

            else:
                if isinstance(attr, obj_model.ManyToManyAttribute):
                    for v in list(secondary_val):
                        secondary_val.remove(v)
                        if v in shared_objs:
                            primary_val.append(v)
                        else:
                            primary_val.append(shared_objs[v])

                elif isinstance(attr, obj_model.OneToManyAttribute):
                    raise ValueError('Unsupported merge operation')

                else:
                    if primary_val and secondary_val and shared_objs.get(secondary_val, None) != primary_val:
                        raise ValueError('Unsupported merge operation')

    def append_obj(self, secondary_obj, shared_objs):
        """ Append secondary object to primary model

        Args:
            secondary_obj (:obj:`obj_model.Model`): object to append to primary model
            shared_objects (:obj:`dict`): dictionary that maps objects in the secondary model to objects in the primary
                model that have the same primary attribute (e.g. id)

        Raises:
            :obj:`ValueError`: if the merge operation is not supported
        """
        for attr_name, attr in secondary_obj.Meta.attributes.items():
            if isinstance(attr, obj_model.RelatedAttribute):
                val = getattr(secondary_obj, attr_name)

                if isinstance(attr, obj_model.ManyToManyAttribute) or \
                        isinstance(attr, obj_model.OneToManyAttribute):
                    for v in list(val):
                        if v in shared_objs:
                            val.remove(v)
                            val.append(shared_objs[v])

                elif isinstance(attr, obj_model.ManyToOneAttribute):
                    if val in shared_objs:
                        secondary_obj.__setattr__(attr_name, shared_objs[val])

                else:
                    raise ValueError('Unsupported merge operation')


    def merge_attributes(self, primary_model, secondary_model, attr):
        """ Merge attributes of two models by merging objects from secondary model into primary model

        Args:
            primary_model (:obj:`Model`): base for merged model
            secondary_model (:obj:`Model`): model to merge into primary model
            attr (:obj:`obj_model.RelatedAttribute`): attribute of models to merge

        Raises:
            :obj:`ValueError`: if merge operation for attribute is not supported
        """
        merge_type = attr.primary_class.Meta.merge
        if merge_type == MergeType.primary:
            continue
        elif merge_type == MergeType.concatenate:
            self.concatenate_attr(primary_model, secondary_model, attr)
        elif merge_type == MergeType.join:
            self.join_attr(primary_model, secondary_model, attr)
        else:
            raise ValueError('Unsupported merge operation: "{}"'.format(merge_type.name)
                             )  # pragma: no cover # unreachable because all values of MergeType are enumerated

    def merge_nested_attributes(self, primary_model, secondary_model, attr):
        """ Merge nested attributes of two models by merging objects from secondary model into primary model

        Args:
            primary_model (:obj:`Model`): base for merged model
            secondary_model (:obj:`Model`): model to merge into primary model
            attr (:obj:`obj_model.RelatedAttribute`): attribute of models to merge

        Raises:
            :obj:`ValueError`: if merge operation for attribute is not supported
        """
        merge_type = attr.primary_class.Meta.merge
        if merge_type == MergeType.primary:
            continue
        elif merge_type == MergeType.concatenate:
            continue
        elif merge_type == MergeType.join:
            self.join_nested_attr(primary_model, secondary_model, attr)
        else:
            raise ValueError('Unsupported merge operation: "{}"'.format(merge_type.name)
                             )  # pragma: no cover # unreachable because all values of MergeType are enumerated

    def concatenate_attr(self, primary_model, secondary_model, attr):
        """ Concatenate the values of the attributes of two models

        Args:
            primary_model (:obj:`Model`): base for merged model
            secondary_model (:obj:`Model`): model to merge into primary model
            attr (:obj:`obj_model.RelatedAttribute`): attribute of models to merge

        Raises:
            :obj:`ValueError`: if the models have objects with the same primary attributes (e.g. ids)
        """
        primary_objs, secondary_objs, shared_objs, different_secondary_objs = self.get_attr_objects(primary_model, secondary_model, attr)

        # verify models have disjoint objects
        if shared_objs:
            raise ValueError('Models {} and {} cannot have {} with the same primary attributes:\n  {}',
                             primary_model.id, secondary_model.id, attr.primary_class.Meta.verbose_name_plural.lower(),
                             '\n  '.join(o.get_primary_attribute() for o in shared_objs))

        # merge objects from secondary model into primary model
        primary_objs.extend(secondary_objs)

    def join_attr(self, primary_model, secondary_model, attr):
        """ Join the values of the attributes of two models

        Args:
            primary_model (:obj:`Model`): base for merged model
            secondary_model (:obj:`Model`): model to merge into primary model
            attr (:obj:`obj_model.RelatedAttribute`): attribute of models to merge

        Raises:
            :obj:`ValueError`: if objects with same prrimary attribute (e.g. id) have different attributes
        """
        primary_objs, secondary_objs, shared_objs, different_secondary_objs = self.get_attr_objects(primary_model, secondary_model, attr)

        # verify that shared objects have same attributes
        for primary_obj, secondary_obj in shared_objs:
            for obj_attr in attr.primary_model.Meta.attributes.values():
                if not isinstance(obj_attr, obj_model.RelatedAttribute):
                    primary_val = getattr(primary_obj, obj_attr.name)
                    secondary_val = getattr(secondary_obj, obj_attr.name)
                    if primary_val != secondary_val:
                        raise ValueError('{} "{}" of {} and {} must have the same value of {}'.format(
                            attr.primary_class.Meta.verbose_name, primary_obj.get_primary_attribute(),
                            primary_model.id, secondary_model.id, obj_attr.name))

        # merge objects from secondary model that aren't shared with primary model
        primary_objs.extend(different_secondary_objs)

    def join_nested_attr(self, primary_model, secondary_model, attr):
        """ Join the values of the attributes of two models

        Args:
            primary_model (:obj:`Model`): base for merged model
            secondary_model (:obj:`Model`): model to merge into primary model
            attr (:obj:`obj_model.RelatedAttribute`): attribute of models to merge

        Raises:
            :obj:`ValueError`: if objects with same prrimary attribute (e.g. id) have different attributes
        """
        _, _, shared_objs, _ = self.get_attr_objects(primary_model, secondary_model, attr)

        # merge objects with same primary attributes
        for primary_obj, secondary_obj in shared_objs:
            for obj_attr_name, obj_attr in itertools.chain(attr.primary_model.Meta.attributes.items(),
                                                           attr.primary_model.Meta.related_attributes.items()):
                if isinstance(obj_attr, obj_model.RelatedAttribute):
                    primary_val = getattr(primary_obj, obj_attr_name)
                    secondary_val = getattr(secondary_obj, obj_attr_name)
                    for v in secondary_val:
                        adf

                    primary_val.extend(secondary_val)

    def get_attr_objects(self, primary_model, secondary_model, attr):
        """ Get the attributes of the primary and secondary models and their intersection and difference

        Args:
            primary_model (:obj:`Model`): base for merged model
            secondary_model (:obj:`Model`): model to merge into primary model
            attr (:obj:`obj_model.RelatedAttribute`): attribute of models to merge

        Returns:
            :obj:`list` of :obj:`obj_model.Model`: value of attribute of primary model
            :obj:`list` of :obj:`obj_model.Model`: value of attribute of secondary model
            :obj:`list` of :obj:`tuple` of :obj:`obj_model.Model`: pairs of objects from the primary and
                secondary models with same primary attribute (e.g. id)
            :obj:`list` of :obj:`obj_model.Model`: objects from the secondary model which have different
                primary attributes (e.g. id) than objects from the primary model
        """
        # get primary and secondary objects
        primary_objs = getattr(primary_model, attr.related_name)
        secondary_objs = getattr(secondary_model, attr.related_name)

        # get shared objects
        primary_vals = [o.get_primary_attribute() for o in primary_objs]
        secondary_vals = [o.get_primary_attribute() for o in secondary_objs]

        shared_vals = set(secondary_vals).intersection(set(secondary_vals))
        shared_objs = []
        for shared_val in shared_vals:
            shared_objs.append((primary_objs[primary_vals.index(shared_val)],
                                secondary_objs[secondary_vals.index(shared_val)]))

        # get differrent objects
        different_secondary_vals = set(secondary_vals).difference(set(primary_vals))
        different_secondary_objs = [secondary_objs[secondary_vals.index(v)] for v in different_secondary_vals]

        # return objects
        return (primary_objs, secondary_objs, shared_objs, different_secondary_objs)

from __future__ import annotations

import glob
import json
import os
import pathlib
import shutil
from enum import StrEnum
from typing import List, Callable, Iterable, TypeVar, Generic

from exception import LeafException, ErrorCodes

T = TypeVar('T', bound='Object')


class ObjectType(StrEnum):
    FILE = 'FILE'
    DIRECTORY = 'DIRECTORY'


class Object:
    def __init__(self, path_or_object_name: str, type: ObjectType, parent: Object | None = None):
        """Create an Object instance
        :param Object parent: The parent filesystem object
        :param str path_or_object_name: Path to the object in the filesystem"""

        if not parent and not os.path.isabs(path_or_object_name):
            raise LeafException(ErrorCodes.MUST_BE_ABS_PATH).with_additional_message(
                f'{path_or_object_name} must be an absolute path')

        self.parent = parent
        self.extension = '.'.join(pathlib.Path(path_or_object_name).suffixes)
        self.name = os.path.basename(path_or_object_name)
        self.type = type

        if parent:
            self.full_path = os.path.join(parent.full_path, self.name)
        else:
            self.full_path = path_or_object_name

        if type == ObjectType.DIRECTORY and not self.physically_exists():
            os.mkdir(self.full_path)
        if type == ObjectType.FILE and not self.physically_exists():
            open(self.full_path, 'a').close()

    def physically_exists(self) -> bool:
        return os.path.exists(self.full_path)

    def _ensure_exists(self):
        if self.physically_exists():
            return
        raise LeafException(ErrorCodes.OBJECT_DOES_NOT_EXIST).with_additional_message(
            f'Object {self} was deleted')

    def __str__(self):
        return f'{self.type}: {self.full_path}'

    def __repr__(self):
        return self.__str__()

    def __hash__(self):
        return self.full_path.__hash__()

    def __eq__(self, other: Object) -> bool:
        return self.full_path == other.full_path

    def __ne__(self, other: Object) -> bool:
        return self.full_path != other.full_path


class FileObject(Object):
    def __init__(self, path_or_object_name: str, parent: Object | None = None):
        super().__init__(path_or_object_name, ObjectType.FILE, parent)

    def read_string(self) -> str:
        self._ensure_exists()

        with open(self.full_path, "r") as f:
            return f.read()

    def read_lines(self) -> List[str]:
        self._ensure_exists()

        with open(self.full_path, "r") as f:
            return f.readlines()

    def read_bytes(self) -> bytes:
        self._ensure_exists()

        with open(self.full_path, "rb") as f:
            return f.read()

    def write_string(self, s: str):
        self._ensure_exists()

        with open(self.full_path, 'w') as f:
            f.write(s)

    def write_bytes(self, b: bytes):
        self._ensure_exists()

        with open(self.full_path, 'wb') as f:
            f.write(b)

    def write_json(self, j):
        self._ensure_exists()

        json.dump(j, open(self.full_path, 'w'))

    def read(self, loader) -> Any:
        self._ensure_exists()

        return loader(open(self.full_path, 'r'))

    def read_json(self) -> T:
        return self.read(json.load)

    def delete(self):
        if self.physically_exists():
            os.remove(self.full_path)


class FileObjects:
    def __init__(self, parent: DirectoryObject):
        self.parent = parent

    def __iter__(self):
        for child in self.parent.where(lambda x: isinstance(x, FileObject)):
            yield child

    def where(self, predicate: Callable[[FileObject], bool], recursive=False) -> Iterable[FileObject]:
        for child in self.parent.children:
            if isinstance(child, FileObject) and predicate(child):
                yield child
            if recursive and isinstance(child, DirectoryObject):
                yield from child.files.where(predicate)

    def first_or_default(self, predicate: Callable[[FileObject], bool], recursive=False) -> FileObject | None:
        """Returns the first or None object in the tree by predicate"""
        for child in self.parent.children:
            if isinstance(child, FileObject) and predicate(child):
                return child
            if recursive and isinstance(child, DirectoryObject):
                found = child.files.first_or_default(predicate, recursive=recursive)
                if found:
                    return found
        return None


class DirectoryObjects:
    def __init__(self, parent: DirectoryObject):
        self.parent = parent

    def __iter__(self):
        for child in self.parent.where(lambda x: isinstance(x, DirectoryObject)):
            yield child

    def where(self, predicate: Callable[[DirectoryObject], bool], recursive=False) -> Iterable[DirectoryObject]:
        for child in self.parent.children:
            if isinstance(child, DirectoryObject) and predicate(child):
                yield child
                if recursive:
                    yield from child.directories.where(predicate)

    def first_or_default(self, predicate: Callable[[DirectoryObject], bool], recursive=False) -> DirectoryObject | None:
        """Returns the first or None object in the tree by predicate"""
        for child in self.parent.children:
            if isinstance(child, DirectoryObject):
                if predicate(child):
                    return child
                if recursive:
                    found = child.directories.first_or_default(predicate, recursive=recursive)
                    if found:
                        return found
        return None


class DirectoryObject(Object):
    def __init__(self, path_or_object_name: str, parent: Object | None = None):
        super().__init__(path_or_object_name, ObjectType.DIRECTORY, parent)

    @property
    def descendants(self) -> Iterable[Object]:
        self._ensure_exists()

        if not self.physically_exists():
            return []

        for path in glob.glob(os.path.join(self.full_path, '*')):
            if os.path.isdir(path):
                d = DirectoryObject(
                    os.path.basename(path),
                    self
                )
                yield d
                yield from d.descendants
            if os.path.isfile(path):
                yield FileObject(
                    os.path.basename(path),
                    self
                )

    @property
    def files(self) -> FileObjects:
        return FileObjects(self)

    @property
    def directories(self) -> DirectoryObjects:
        return DirectoryObjects(self)

    @property
    def children(self) -> Iterable[Object]:
        self._ensure_exists()

        if not self.physically_exists():
            return []

        for path in glob.glob(os.path.join(self.full_path, '*')):
            if os.path.isdir(path):
                yield DirectoryObject(
                    os.path.basename(path),
                    self
                )
            if os.path.isfile(path):
                yield FileObject(
                    os.path.basename(path),
                    self
                )

    def add_file(self, name: str) -> FileObject:
        o = FileObject(
            name,
            self
        )

        return o

    def add_directory(self, name: str) -> DirectoryObject:
        o = DirectoryObject(
            name,
            self
        )

        return o

    def delete(self):
        if self.physically_exists():
            shutil.rmtree(self.full_path)

    def where(self, predicate: Callable[[Object], bool], recursive=False) -> Iterable[Object]:
        self._ensure_exists()

        for child in self.children:
            if predicate(child):
                yield child
            if recursive and isinstance(child, DirectoryObject):
                yield from child.where(predicate)

    def exists(self, predicate: Callable[[Object], bool], recursive=False) -> bool:
        """Check if the object exists in the tree."""
        self._ensure_exists()

        for child in self.children:
            if predicate(child):
                return True
            if recursive and isinstance(child, DirectoryObject) and child.exists(predicate, recursive=recursive):
                return True

        return False

    def first_or_default(self, predicate: Callable[[Object], bool], recursive=False) -> Object | None:
        """Returns the first or None object in the tree by predicate"""
        self._ensure_exists()

        for child in self.children:
            if predicate(child):
                return child
            if recursive and isinstance(child, DirectoryObject):
                found = child.first_or_default(predicate, recursive=recursive)
                if found:
                    return found
        return None

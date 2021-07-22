



// export type MongooseFilterQuery<T> = {
//     [P in keyof T]?: P extends '_id'
//         ? [Extract<T[P], mongodb.ObjectId>] extends [never]
//             ? mongodb.Condition<T[P] | string | (T[P] | string)[]>
//             : mongodb.Condition<T[P] | string | { _id: mongodb.ObjectId }>
//         : [Extract<T[P], mongodb.ObjectId>] extends [never]
//             ? mongodb.Condition<T[P] | string | (T[P] | string)[]>
//             : mongodb.Condition<T[P] | string | (T[P] | string)[]>
// } &
//     mongodb.RootQuerySelector<T>;

//export type CBDocument<T> = mongoose.Document & T;

// Remove key from type
export type TypeWithoutKey<Type, Key> = Pick<Type, Exclude<keyof Type, Key>>;










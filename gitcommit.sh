echo "========================================="
read -p "please input commit message:" message
git add .
git commit -m "$message"
echo "DONE !!!"
echo "========================================="
